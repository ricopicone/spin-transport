from spin_transport_01 import simulate
import numpy as np
from matplotlib import pyplot as plt
from bc_generator import dirichlet_bc
import matplotlib.backends.backend_pdf

with matplotlib.backends.backend_pdf.PdfPages("oscillation.pdf") as pdf:

	num_cells = np.logspace(2.5, 3, 2)

	plt.figure()

	for cells in num_cells:
		try:
			result = simulate(
				Bloch = False,
				Pulse = False,
				DirichletBCleft = dirichlet_bc(-7.5),
				DirichletBCright = dirichlet_bc(7.5),
				n = int(cells))
		except RuntimeError as e:
			print('n = {}'.format(int(cells)))
			raise e
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
		try:
			result = simulate(
				Bloch = False,
				Pulse = False,
				DirichletBCleft = dirichlet_bc(-7.5),
				DirichletBCright = dirichlet_bc(7.5),
				num_steps = int(steps),
				n = 500)
		except RuntimeError as e:
			print('num_steps = {}'.format(int(steps)))
			raise e
		rho1 = result['rho1'][-1,:]
		x = result['x']
		plt.plot(x, rho1, label = 'n = {}'.format(int(steps)))

	plt.legend()
	plt.title('Steady State Solution of $\\rho_1$ with Variable Number of Timesteps')
	plt.xlabel('$\\bar{r}$')
	plt.ylabel('$\\rho_1$')

	pdf.savefig(plt.gcf())

	pdf.close()
