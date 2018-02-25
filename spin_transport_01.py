from fenics import *
import numpy as np
import physical_constants_etc as constants
import argparse

# Parameters
gamma2 = 1 # gyromagnetic ratio for nuclei
gamma3 = 1 # gyromagnetic ratio for electrons
Gamma2 = 1 # transport coefficients for nuclear polarization
Gamma3 = 1 # transport coefficients for electron polarization
Delta2 = 1 # constant spin densities in solid medium for nuclei
Delta3 = 1 # onstant spin densities in solid medium for electrons
Bd = 1 # maximum dipole-dipole magnetic field (Is this a constant?)
B1p = 1 # (What is this and is it a constant?)
B1e = 1 # (What is this and is it a constant?)
T1p = 1 # (What is this and is it a constant?)
T1e = 1 # (What is this and is it a constant?)
T2p = 1 # (What is this and is it a constant?)
T2e = 1 # (What is this and is it a constant?)
rho20 = 1 # (What is this and is it a constant?)
rho30 = 1 # (What is this and is it a constant?)
mu2 = 1 # (What is this and is it a constant?)
mu3 = 1 # (What is this and is it a constant?)
kB = 100 # (What is this and is it a constant?)
Temp = 1 # Temperature used in SS solution

# Defining an invese hyperbolic tangent because fenics doesn't have one
def atanh(x): # Only good for abs(x) < 1!
	return ln((1 + x) / (1 - x)) / 2

def simulate(
	T = 3.0, # final time
	num_steps = 100, # number of time steps
	L = 1, # The lendth of the mesh
	n = 50, # The number of cells in the mesh
	DirichletBCleft = None, # The left Dirichlet BC values
	DirichletBCright = None, # The right Dirichlet BC values
	NeumannBC = None, # The Nuemann BC values
	BSteadyState = 'x[0]' # The steady state Magnetic field expression
	):

	dt = T / num_steps # time step size

	# Initialize data arrays
	x = np.linspace(0, L, n + 1)
	t = np.arange(dt, T + dt, dt)
	_rho1 = np.ndarray([num_steps, n + 1])
	_rho2 = np.ndarray([num_steps, n + 1])
	_rho3 = np.ndarray([num_steps, n + 1])

	# Define the mesh
	mesh = IntervalMesh(n, 0, L) # A 1d mesh with n cells from 0 to L

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
	B = Expression('x[0]', degree = 1)

	# Use an expression to get r
	r = Expression('x[0]', degree = 1)

	# Simulation Constants
	_dt = Constant(dt)
	dr = dx
	gamma = Constant(gamma3 / gamma2)
	Gamma = Constant(Gamma3 / Gamma2)
	Delta = Constant(Delta3 / Delta2)

	# Initial Conditions
	rho_ic = Expression(
		(
			'tanh(((1 + gamma * Delta) / (gamma * (1 + Delta))) * (mu3 * Bd / (kB * Temp)))',
			'tanh(mu2 * ({}) / (kB * Temp))'.format(BSteadyState),
			'tanh(mu3 * ({}) / (kB * Temp))'.format(BSteadyState)),
		gamma = gamma,
		Delta = Delta,
		mu2 = mu2,
		mu3 = mu3,
		Bd = Bd,
		Temp = Temp,
		kB = kB,
		degree = 1)

	# Set initial conditions
	rho.assign(rho_ic)
	rho_n.assign(rho_ic)

	# Functions
	c = B * (1 + Delta) / (1 + gamma * Delta)
	w1p = -gamma * Constant(B1p) # Is this gamma bar?
	w1e = -gamma * Constant(B1e) # Is this gamma bar?
	delta = gamma * Constant(Bd) * r # Is this really gamma bar?
	tau_p = 1 / (1 / T1p + (T2p * w1p**2 / (1 + T2p**2 * delta**2))) # I don't think these should be the same
	tau_e = 1 / (1 / T1e + (T2e * w1e**2 / (1 + T2e**2 * delta**2))) # I don't think these should be the same

	# Define variational problem
	F = ((rho1 - rho1_n) * v1 / _dt
	  + (1 + Gamma) * grad(rho1)[0] * grad(v1)[0]
	  - (c**2 / (1 + Delta)) * ((1 - rho2**2) + Gamma * Delta * gamma**2 * (1 - rho3**2)) * v1 * atanh(rho1)
	  - (c / (1 + Delta) * (grad(rho2)[0] - Gamma * Delta * gamma * grad(rho3)[0]) * v1)) * dr \
	  + ((rho2 - rho2_n) * v2 / _dt
	  + grad(rho2)[0] * grad(v2)[0]
	  + (-c * ((1 - rho2**2) / (1 - rho1**2)) * grad(rho1)[0]
	  + 2 * c * rho2 * atanh(rho1) * grad(rho2)[0] + rho2 / tau_p) * v2) * dr \
	  + ((rho3 - rho3_n) * v3 / _dt
	  + Gamma * grad(rho3)[0] * grad(v3)[0]
	  + (-c * ((1 - rho3**2) / (1 - rho1**2)) * grad(rho1)[0]
	  + 2 * c * rho3 * atanh(rho1) * grad(rho3)[0]
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

	# Time-stepping
	_t = 0
	for n in range(num_steps):

		# Update current time
		_t += dt

		# Solve variational problem for time step
		solve(F == 0, rho, bcs = bc)

		# Add the data
		_rho = rho.vector().get_local().reshape([x.shape[0], len(rho.split())])
		_rho1[n,:] = _rho[:,0]
		_rho2[n,:] = _rho[:,1]
		_rho3[n,:] = _rho[:,2]

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
	args = parser.parse_args()

	data = simulate()

	np.savez(args.s, **data)
