from spin_transport_01 import simulate
import numpy as np
from matplotlib import pyplot as plt
from bc_generator import dirichlet_bc
import matplotlib.backends.backend_pdf

with matplotlib.backends.backend_pdf.PdfPages("oscillation.pdf") as pdf:

	num_cells = np.logspace(2.5, 4.5, 5)

	plt.figure()

	for cells in num_cells:
		result = simulate(
			Bloch = False,
			Pulse = False,
			DirichletBCleft = dirichlet_bc(-7.5),
			DirichletBCright = dirichlet_bc(7.5),
			n = int(cells))
		rho1 = result['rho1'][-1,:]
		x = result['x']
		plt.plot(x, rho1, label = 'n = {}'.format(int(cells)))

	plt.legend()
	plt.title('Steady State Solution of $\\rho_1$ with Variable Number of Cells')
	plt.xlabel('$\\bar{r}$')
	plt.ylabel('$\\rho_1$')

	pdf.savefig(plt.gcf())

	num_steps = np.logspace(1, 2.5, 4)

	plt.figure()

	for steps in num_steps:
		result = simulate(
			Bloch = False,
			Pulse = False,
			DirichletBCleft = dirichlet_bc(-7.5),
			DirichletBCright = dirichlet_bc(7.5),
			num_steps = int(steps),
			n = 500)
		rho1 = result['rho1'][-1,:]
		x = result['x']
		plt.plot(x, rho1, label = 'n = {}'.format(int(steps)))

	plt.legend()
	plt.title('Steady State Solution of $\\rho_1$ with Variable Number of Timesteps')
	plt.xlabel('$\\bar{r}$')
	plt.ylabel('$\\rho_1$')

	pdf.savefig(plt.gcf())

	lengths = np.linspace(7, 8, 200)
	mag = []

	plt.figure()

	for l in lengths:
		try:
			result = simulate(
				Bloch = False,
				Pulse = False,
				DirichletBCleft = dirichlet_bc(-l / 2.),
				DirichletBCright = dirichlet_bc(l / 2.),
				n = 700,
				L = l)
		except RuntimeError as e:
			print('length = {}'.format(l))
			raise e
		rho1 = result['rho1'][-1,:]
		mag.append(np.max(rho1) - np.min(rho1))

	plt.semilogy(lengths, mag)

	plt.title('$\\rho_1$ Oscillation Magnitude vs Mesh Length')
	plt.xlabel('$L$')
	plt.ylabel('Osicllation Magnitude')

	pdf.savefig(plt.gcf())
