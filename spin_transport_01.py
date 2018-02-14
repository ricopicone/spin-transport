from fenics import *

T = 2.0            # final time
num_steps = 1000    # number of time steps
dt = T / num_steps # time step size

dbcl = [0.5, 0.2, 0.3] # The value of the left Dirichlet BC
dbcr = [0.5, 0.1, 0.4] # The value of the right Dirichlet BC

L = 1	# The length of the mesh
n = 200	# The number of cells in the mesh

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

# Create VTK files for visualization output
vtkfile_rho1 = File('spin_transport_soln/rho_1.pvd')
vtkfile_rho2 = File('spin_transport_soln/rho_2.pvd')
vtkfile_rho3 = File('spin_transport_soln/rho_3.pvd')

# Create progress bar
progress = Progress('Time-stepping')
set_log_level(PROGRESS)

# Time-stepping
t = 0
for n in range(num_steps):

	# Update current time
	t += dt

	# Solve variational problem for time step
	solve(F == 0, rho, bcs = dbc)

	# Save solution to file (VTK)
	_rho1, _rho2, _rho3 = rho.split()
	vtkfile_rho1 << (_rho1, t)
	vtkfile_rho2 << (_rho2, t)
	vtkfile_rho3 << (_rho3, t)

	# Update previous solution
	rho_n2.assign(rho_n1)
	rho_n1.assign(rho)

	# Update progress bar
	progress.update(t / T)
