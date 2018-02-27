from fenics import *
import numpy as np
from physical_constants_etc import *
from experimental_constants import *
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
	L = 2e-5, # The lendth of the mesh
	n = 50, # The number of cells in the mesh
	DirichletBCleft = None, # The left Dirichlet BC values
	DirichletBCright = None, # The right Dirichlet BC values
	NeumannBC = None, # The Nuemann BC values
	quiet = False # Suppress log output
	):

	dt = T / num_steps # time step size

	# Initialize data arrays
	x = np.linspace(L/2, -L/2, n + 1)
	t = np.arange(0, T + dt, dt)
	_rho1 = np.ndarray([num_steps + 1, n + 1])
	_rho2 = np.ndarray([num_steps + 1, n + 1])
	_rho3 = np.ndarray([num_steps + 1, n + 1])

	# Define the mesh
	mesh = IntervalMesh(n, -L/2, L/2) # A 1d mesh with n cells from 0 to L

	# Define function space for system of concentrations
	V = VectorFunctionSpace(mesh, 'CG', 1, dim = 3)

	# Setup the Dirichlet BC
	left = lambda x, on_boundary: near(x[0], 0) and on_boundary
	right = lambda x, on_boundary: near(x[0], L) and on_boundary
	bc = []
	try:
		bc.append(DirichletBC(V, Constant(DirichletBCleft), left))
	except:
		pass
	try:
		bc.append(DirichletBC(V, Constant(DirichletBCright), right))
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

	# Set up Magnetic Field Function
	B = Expression('grad * x[0]', grad = Grad, degree = 1)

	# Use an expression to get r
	r = Expression('x[0]', degree = 1)

	# Initial Conditions
	rho_ic = Expression(
		(
			'tanh(((1 + gamma * Delta) / (gamma * (1 + Delta))) * (mu3 * Bd / (kB * Temp)))',
			'tanh(mu2 * (grad * x[0]) / (kB * Temp))',
			'tanh(mu3 * (grad * x[0]) / (kB * Temp))'),
		gamma = γb,
		Delta = δb,
		mu2 = μe,
		mu3 = μp,
		Bd = Bd,
		Temp = temp,
		kB = kB,
		grad = Grad,
		degree = 1)
	_, rho20, rho30 = split(rho_ic)

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

	# Define variational problem
	F = ((rho1 - rho1_n) * v1 / _dt
	  + (1 + Γb) * grad(rho1)[0] * grad(v1)[0]
	  - (cb**2 / (1 + δb)) * ((1 - rho2**2) + Γb * δb * γb**2 * (1 - rho3**2)) * v1 * atanh(rho1)
	  - (cb / (1 + δb) * (grad(rho2)[0] - Γb * δb * γb * grad(rho3)[0]) * v1)) * dr \
	  + ((rho2 - rho2_n) * v2 / _dt
	  + grad(rho2)[0] * grad(v2)[0]
	  + (-cb * ((1 - rho2**2) / (1 - rho1**2)) * grad(rho1)[0]
	  + 2 * cb * rho2 * atanh(rho1) * grad(rho2)[0] + rho2 / tau_p) * v2) * dr \
	  + ((rho3 - rho3_n) * v3 / _dt
	  + Γb * grad(rho3)[0] * grad(v3)[0]
	  + (-cb * ((1 - rho3**2) / (1 - rho1**2)) * grad(rho1)[0]
	  + 2 * cb * rho3 * atanh(rho1) * grad(rho3)[0]
	  + rho3 / tau_e) * v3) * dr \
	  - (rho20 * v2 / T1p + rho30 * v3 / T1e) * dr

	# Add the Neumann BC
	if NeumannBC is not None:
		G1 = Constant(NeumannBC[0])
		G2 = Constant(NeumannBC[1])
		G3 = Constant(NeumannBC[2])
		F += (v1 * G1 + v2 * G2 + v3 * G3) * ds

	# Create progress bar
	progress = Progress('Time-stepping')
	set_log_level(PROGRESS)
	set_log_active(not quiet)

	# Add the data
	_rho = rho.vector().get_local().reshape([x.shape[0], len(rho.split())])
	_rho1[0,:] = _rho[:,0]
	_rho2[0,:] = _rho[:,1]
	_rho3[0,:] = _rho[:,2]

	# Time-stepping
	_t = 0
	for i in range(1, num_steps + 1):

		# Update current time
		_t += dt

		# Solve variational problem for time step
		solve(F == 0, rho, bcs = bc)

		# Add the data
		_rho = rho.vector().get_local().reshape([x.shape[0], len(rho.split())])
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
		default = 2e-5)
	parser.add_argument('-n',
		type = int,
		help = 'The number of cells to use in the mesh.',
		default = 50)
	args = parser.parse_args()

	data = simulate(
		T = args.t,
		num_steps = args.T,
		L = args.L,
		n = args.n)

	np.savez(args.s, **data)
