from fenics import *

T = 5.0            # final time
num_steps = 500    # number of time steps
dt = T / num_steps # time step size

L = 1	# The length of the mesh
n = 20	# The number of cells in the mesh

mesh = IntervalMesh(n, 0, L) # A 1d mesh with n cells from 0 to L

# Define function space for system of concentrations
#element = VectorElement('P', triangle, 1, dim = 3) # Vector of three linear lagrangian finite elements for the three variable we are trying to solve for
#V = FunctionSpace(mesh, element) # Creating the function space
V = VectorFunctionSpace(mesh, 'CG', 1, dim = 3)

# Define test functions
v_1, v_2, v_3 = TestFunctions(V)

rho = Function(V)
rho_n1 = Function(V)
rho_n2 = Function(V)

# Split system functions to access components
rho1, rho2, rho3 = split(rho)
rho1_n1, rho2_n1, rho3_n1 = split(rho_n1)
rho1_n2, rho2_n2, rho3_n2 = split(rho_n2)

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
	solve(F == 0, rho)

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
