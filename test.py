import unittest
from spin_transport_01 import *
import numpy as np	
	
class TestSSSolution(unittest.TestCase):

	def sim_test(self, solution, tol, **kwargs):
		data = simulate(**kwargs)
		result = True
		for key in solution:
			result = result and np.all(np.abs(solution[key](data['x']) - data[key]) < tol)
		if not result:
			print('plot error')
		self.assertTrue(result)

	def test_SS(self):
		self.sim_test(
			{
				'rho1': lambda x: 0.1,
				'rho2': lambda x: np.tanh(x),
				'rho3': lambda x: np.tanh(x)
			},
			0.01)

if __name__ == '__main__':
	unittest.main()
