import unittest
from spin_transport_01 import *
import numpy as np

class TestSSSolution(unittest.TestCase):

	def sim_test(self, solution, tol, **kwargs):
		kwargs['quiet'] = True
		data = simulate(**kwargs)
		for key in solution:
			self.assertTrue(
				np.all(np.abs(
					solution[key](data['x']) - data[key][-1,:]) < tol),
				'{} failed.'.format(key))

	def test_SS(self):
		self.sim_test(
			{
				'rho1': lambda x: np.tanh(((1 + γb * δb) / (γb * (1 + δb))) * (μp * Bd / (kB * temp))),
				'rho2': lambda x: np.tanh(μe * (Grad * x) / (kB * temp)),
				'rho3': lambda x: np.tanh(μp * (Grad * x) / (kB * temp))
			},
			0.005)

if __name__ == '__main__':
	unittest.main()
