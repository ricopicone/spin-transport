from fenics import *
import numpy as np
import physical_constants_etc as constants
import argparse

parser = argparse.ArgumentParser(
	description = 'Spin Transport Simulation',
	add_help = False)
parser.add_argument('-s',
	type = str,
	help = 'The file to save the solution to.',
	default = 'spin_transport_soln/soln.npz')
args = parser.parse_args()

T = 2.0            # final time
num_steps = 500    # number of time steps
dt = T / num_steps # time step size

dbcl = [0.5, 0.2, 0.3] # The value of the left Dirichlet BC
dbcr = [0.5, 0.1, 0.4] # The value of the right Dirichlet BC

L = 1	# The length of the mesh
n = 100	# The number of cells in the mesh

# Initialize data arrays
x = np.linspace(0, L, n + 1)
t = np.arange(dt, T + dt, dt)
_rho1 = np.ndarray([num_steps, n + 1])
_rho2 = np.ndarray([num_steps, n + 1])
_rho3 = np.ndarray([num_steps, n + 1])

# Initial Conditions
rho1_ic = Expression(('0.5 * sin(pi * x[0])', '0.15 * sin(pi * x[0])', '0.35 * sin(pi * x[0])'), degree = 1)

mesh = IntervalMesh(n, 0, L) # A 1d mesh with n cells from 0 to L

# Define function space for system of concentrations
V = VectorFunctionSpace(mesh, 'CG', 1, dim = 3)

# Define the boundary conditions
def left(x, on_boundary):
	return near(x[0], 0) and on_boundary
def right(x, on_boundary):
	return near(x[0], L) and on_boundary

dbc = [
	DirichletBC(V, Constant(dbcl), left),
	DirichletBC(V, Constant(dbcr), right)
]

# Define test functions
v_1, v_2, v_3 = TestFunctions(V)

rho = Function(V)
rho_n1 = Function(V)
rho_n2 = Function(V)

# Split system functions to access components
rho1, rho2, rho3 = split(rho)
rho1_n1, rho2_n1, rho3_n1 = split(rho_n1)
rho1_n2, rho2_n2, rho3_n2 = split(rho_n2)

# Set initial conditions
rho.assign(rho1_ic)
rho_n1.assign(rho1_ic)
rho_n2.assign(rho1_ic)

# Define variational problem
_dt = Constant(dt)
dr = dx
F = ((rho1 - 2 * rho1_n1 + rho1_n2) / _dt**2) * v_1 * dr \
  + ((rho1 - rho1_n1) / _dt) * nabla_grad(v_1)[0] * dr \
  + ((rho2 - 2 * rho2_n1 + rho2_n2) / _dt**2) * v_2 * dr \
  + ((rho2 - rho2_n1) / _dt) * nabla_grad(v_2)[0] * dr \
  + ((rho3 - 2 * rho3_n1 + rho3_n2) / _dt**2) * v_3 * dr \
  + ((rho3 - rho3_n1) / _dt) * nabla_grad(v_3)[0] * dr

# Create progress bar
progress = Progress('Time-stepping')
set_log_level(PROGRESS)

# Time-stepping
_t = 0
for n in range(num_steps):

	# Update current time
	_t += dt

	# Solve variational problem for time step
	solve(F == 0, rho, bcs = dbc)

	# Add the data
	_rho = rho.vector().get_local().reshape([x.shape[0], len(rho.split())])
	_rho1[n,:] = _rho[:,0]
	_rho2[n,:] = _rho[:,1]
	_rho3[n,:] = _rho[:,2]

	# Update previous solution
	rho_n2.assign(rho_n1)
	rho_n1.assign(rho)

	# Update progress bar
	progress.update(_t / T)

np.savez(args.s, rho1 = _rho1, rho2 = _rho2, rho3 = _rho3, t = t, x = x)
