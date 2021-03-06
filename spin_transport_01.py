from fenics import *
import numpy as np
from physical_constants_etc import *
from experimental_constants import *
from bc_generator import *
import argparse

# Other constants
B1e = 0.05
B1p = 1

# Defining an invese hyperbolic tangent because fenics doesn't have one
def atanh(x): # Only good for abs(x) < 1!
	return ln((1 + x) / (1 - x)) / 2

def simulate(
	T = 0.1, # final time
	num_steps = 50, # number of time steps
	L = 15, # The length of the mesh
	n = 50, # The number of cells in the mesh
	DirichletBCleft = None, # The left Dirichlet BC values
	DirichletBCright = None, # The right Dirichlet BC values
	NeumannBC = None, # The Nuemann BC values
	Seperation = True, # Enable the Seperation
	Bloch = True, # Enable the Bloch Dynamics
	Diffusion = True, # Enable Diffusion
	Pulse = True, # Enable Bloch Pulse Dynamics
	ic = None, # The system initial conditions
	fail = True, # Whether to throw an error when the solver fails
	quiet = False # Suppress log output
	):

	ds = globals()['ds']

	dt = T / num_steps # time step size

	# Initialize data arrays
	x = np.linspace(L/2, -L/2, n + 1)
	t = np.arange(0, T + dt, dt)
	_rho1 = np.ndarray([num_steps + 1, n + 1])
	_rho2 = np.ndarray([num_steps + 1, n + 1])
	_rho3 = np.ndarray([num_steps + 1, n + 1])
	_dx = np.gradient(x)[0]

	# Define the mesh
	mesh = IntervalMesh(n, -L/2, L/2) # A 1d mesh with n cells from 0 to L

	# Define function space for system of concentrations
	V = VectorFunctionSpace(mesh, 'CG', 1, dim = 3)

	# Setup the Dirichlet BC
	left = lambda x, on_boundary: near(x[0], -L / 2.) and on_boundary
	right = lambda x, on_boundary: near(x[0], L / 2.) and on_boundary
	bcs = []
	try:
		bcs.append(DirichletBC(V, Constant(DirichletBCleft), left))
	except:
		pass
	try:
		bcs.append(DirichletBC(V, Constant(DirichletBCright), right))
	except:
		pass

	# Define test functions
	v1, v2, v3 = TestFunctions(V)

	# Define trial functions
	rho = Function(V)
	rho_n = Function(V)

	# Split system functions to access components
	rho1, rho2, rho3 = split(rho)
	rho1_n, rho2_n, rho3_n = split(rho_n)

	# Set up Magnetic Field Function -- tilde version!
	Bt = Expression('b0 + bd * x[0]', b0 = B0, bd = Bd, degree = 1)
	# ... because Bt(rb) = B(r) = B0 + g*r = B0 + g*(rb/(g/Bd)) = B0 + Bd*rb
	# ... and, as I understand, rb = x[0]

	# Use an expression to get r
	r = Expression('x[0]', degree = 1)

	if ic is None:
		ic = (
			'-tanh(((1 + gamma * Delta) / (gamma * (1 + Delta))) * (mu3 * Bd / (kB * Temp)))',
			'tanh(mu2 * Bt / (kB * Temp))',
			'tanh(mu3 * Bt / (kB * Temp))')

	# Initial Conditions
	rho_ic = Expression(ic,
		gamma = γb,
		Delta = δb,
		mu2 = μp,
		mu3 = μe,
		Bd = Bd,
		Temp = temp,
		kB = kB,
		grad = Grad,
		Bt = Bt,
		degree = 1)
	#_, rho20, rho30 = split(rho_ic)
	rho20 = rho_ic[1]
	rho30 = rho_ic[2]

	# Set initial conditions
	rho.assign(rho_ic)
	rho_n.assign(rho_ic)

	# Functions
	w1p = -Γb * B1p
	w1e = -Γb * B1e
	delta = Γb * Bd * r
	tau_p = 1 / (1 / T1p + (T2p * w1p**2 / (1 + T2p**2 * delta**2)))
	tau_e = 1 / (1 / T1e + (T2e * w1e**2 / (1 + T2e**2 * delta**2)))

	# Simulation Constants
	_dt = Constant(dt)
	dr = dx

	# Define the Base Variational Problem
	F = ((rho1 - rho1_n) * v1 / _dt
	  + (rho2 - rho2_n) * v2 / _dt
	  + (rho3 - rho3_n) * v3 / _dt) * dr

	# Add the Seperation
	if Seperation:
		print('Seperation Enabled')
		F += -((cb**2 / (1 + δb)) * ((1 - rho2**2) + Γb * δb * γb**2 * (1 - rho3**2)) * v1 * atanh(rho1)
		  - (cb / (1 + δb) * (grad(rho2)[0] - Γb * δb * γb * grad(rho3)[0]) * v1)) * dr \
		  + ((-cb * ((1 - rho2**2) / (1 - rho1**2)) * grad(rho1)[0]
		  + 2 * cb * rho2 * atanh(rho1) * grad(rho2)[0]) * v2) * dr \
		  + ((-cb * ((1 - rho3**2) / (1 - rho1**2)) * grad(rho1)[0]
		  + 2 * cb * rho3 * atanh(rho1) * grad(rho3)[0]) * v3) * dr
		
	# Add the Diffusion
	if Diffusion:
		print('Diffusion Enabled')
		F += (1 + Γb) * grad(rho1)[0] * grad(v1)[0] * dr \
		  + grad(rho2)[0] * grad(v2)[0] * dr \
		  + Γb * grad(rho3)[0] * grad(v3)[0] * dr

	# Add the Bloch Relaxation Dynamics
	if Bloch:
		print('Bloch Relaxation Enabled')
		F += -(rho20 * v2 / T1p + rho30 * v3 / T1e) * dr

	# Add the Bloch Pulse Dynamics
	if Pulse:
		print('Bloch Pulse Dynamics Enabled')
		F += (rho2 / tau_p) * v2 * dr \
		  + (rho3 / tau_e) * v3 * dr

	class LeftEdge(SubDomain):
		def inside(self, x, on_boundary):
			return on_boundary and near(x[0], -L/2.)

	class RightEdge(SubDomain):
		def inside(self, x, on_boundary):
			return on_boundary and near(x[0], L/2.)

	boundaryLeft = LeftEdge()
	boundaryRight = RightEdge()
	
	boundaries = FacetFunction("size_t", mesh)
	boundaryLeft.mark(boundaries,1)
	boundaryRight.mark(boundaries,2)
	ds = ds(subdomain_data = boundaries)

	# Add the Neumann BC
	if NeumannBC is not None:
		G1 = Constant(NeumannBC[0])
		G2 = Constant(NeumannBC[1])
		G3 = Constant(NeumannBC[2])

		F += (v1 * G1 + v2 * G2 + v3 * G3) * ds(1)
		F += (v1 * G1 - v2 * G2 - v3 * G3) * ds(2)

	# Create progress bar
	progress = Progress('Time-stepping')
	set_log_level(PROGRESS)
	set_log_active(not quiet)

	# Add the data
	_rho = rho.vector().get_local().reshape([x.shape[0], len(rho.split())])
	_rho1[0,:] = _rho[:,0]
	_rho2[0,:] = _rho[:,1]
	_rho3[0,:] = _rho[:,2]

	# Define the solver
	J = derivative(F, rho)

	class NLP(NonlinearProblem):
		def F(self, b, x):
			assemble(F, tensor=b)
			for bc in bcs:
				bc.apply(b, x)
		def J(self, A, x):
			assemble(J, tensor=A)
			for bc in bcs:
				bc.apply(A)

	solver = PETScSNESSolver()
	PETScOptions.set("ksp_type", "preonly")
	PETScOptions.set("pc_type", "lu")
	PETScOptions.set("pc_factor_mat_solver_type", "mumps")
	PETScOptions.set("snes_atol", 1e-8)
	solver.set_from_options()

	# Time-stepping
	_t = 0
	for i in range(1, num_steps + 1):

		# Update current time
		_t += dt

		j = 0
		while True:
			# Solve variational problem for time step
			try:
				solver.solve(NLP(), rho.vector())
			except RuntimeError as e:
				if fail:
					raise e
				else:
					print(e)
					print('At step {}: t = {}s'.format(i, _t))
					break

			# Add the data
			_rho = rho.vector().get_local().reshape([x.shape[0], len(rho.split())])

			# If using Neumann BC iteratively fix the BC
			if NeumannBC is not None:
				j += 1
				d_rho = np.gradient(_rho, _dx, axis = 0)
				if j >= 10 or np.max(np.abs(np.array(NeumannBC) - d_rho[0,:])) < 1e-10:
					break
				print('Iteration {}'.format(j))
				print('Gradient rho1: {} rho2: {} rho3: {}'.format(d_rho[0,0], d_rho[0,1], d_rho[0,2]))
				print('BC rho1: {} rho2: {} rho3: {}'.format(*NeumannBC))
				print('Set BC rho1: {} rho2: {} rho3: {}'.format(G1.values()[0], G2.values()[0], G3.values()[0]))
				G1.assign(G1.values()[0] * NeumannBC[0] / d_rho[0,0])
				G2.assign(G2.values()[0] * NeumannBC[1] / d_rho[0,1])
				G3.assign(G3.values()[0] * NeumannBC[2] / d_rho[0,2])
				print('New set BC rho1: {} rho2: {} rho3: {}'.format(G1.values()[0], G2.values()[0], G3.values()[0]))
			else:
				break

		_rho1[i,:] = _rho[:,0]
		_rho2[i,:] = _rho[:,1]
		_rho3[i,:] = _rho[:,2]

		# Update previous solution
		rho_n.assign(rho)

		# Update progress bar
		progress.update(_t / T)

	return {
		'rho1': _rho1,
		'rho2': _rho2,
		'rho3': _rho3,
		't': t,
		'x': x}

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description = 'Spin Transport Simulation')
	parser.add_argument('-s',
		type = str,
		help = 'The file to save the solution to.',
		default = 'spin_transport_soln/soln.npz')
	parser.add_argument('-t',
		type = float,
		help = 'The duration of time to simulate.',
		default = 0.1)
	parser.add_argument('-T',
		type = int,
		help = 'The number of timesteps to take.',
		default = 50)
	parser.add_argument('-L',
		type = float,
		help = 'The length of the mesh.',
		default = 15)
	parser.add_argument('-n',
		type = int,
		help = 'The number of cells to use in the mesh.',
		default = 50)
	parser.add_argument('-S',
		action = 'store_false',
		help = 'Disable Seperation.')
	parser.add_argument('-b',
		action = 'store_false',
		help = 'Disable Bloch Relaxation Dynamics.')
	parser.add_argument('-p',
		action = 'store_false',
		help = 'Disable Bloch Pulse Dynamics.')
	parser.add_argument('-d',
		action = 'store_false',
		help = 'Disable Diffusion')
	parser.add_argument('-f',
		action = 'store_false',
		help = 'Don\'t raise exception on solver failure.')
	parser.add_argument('--NeumannBC',
		type = float,
		nargs = 3,
		help = 'The Neumann Boundary Conditions to apply.',
		metavar = ('grad(rho1)', 'grad(rho2)', 'grad(rho3)'),
		default = None),
	parser.add_argument('--DirichletBCleft',
		type = float,
		nargs = 3,
		metavar = ('rho1', 'rho2', 'rho3'),
		help = 'The Dirichlet Boundary Conditions to apply on the left side.',
		default = None)
	parser.add_argument('--DirichletBCright',
		type = float,
		nargs = 3,
		metavar = ('rho1', 'rho2', 'rho3'),
		help = 'The Dirichlet Boundary Conditions to apply on the right side.',
		default = None)
	parser.add_argument('-i',
		type = str,
		nargs = 3,
		metavar = ('rho1_ic', 'rho2_ic', 'rho3_ic'),
		help = 'The Initial Conditions to use for starting the simulation.',
		default = None)
	parser.add_argument('--DirichletSSBC',
		action = 'store_true',
		help = 'Set Dirichlet Boundary Conditions to the Steady State Initial Conditions.')
	parser.add_argument('--NeumannSSBC',
		action = 'store_true',
		help = 'Set Neumann Boundary Conditions to the Steady State Initial Conditions.')
	args = parser.parse_args()

	dirichlet = [
		args.DirichletBCleft,
		args.DirichletBCright]
	neumann = args.NeumannBC

	if args.DirichletSSBC:
		dirichlet = [
			dirichlet_bc(-args.L / 2.),
			dirichlet_bc(args.L / 2.)]

	if args.NeumannSSBC:
		neumann = neumann_bc(args.L / 2.)
	
	data = simulate(
		T = args.t,
		num_steps = args.T,
		L = args.L,
		n = args.n,
		NeumannBC = neumann,
		DirichletBCleft = dirichlet[0],
		DirichletBCright = dirichlet[1],
		Seperation = args.S,
		Bloch = args.b,
		Diffusion = args.d,
		Pulse = args.p,
		ic = args.i,
		fail = args.f)

	np.savez(args.s, **data)
