import unittest
from spin_transport_01 import *
import numpy as np

class TestSSBaseSys(unittest.TestCase):

	@classmethod
	def setUpClass(self):
		self.ic = [0.1, 0.2, 0.3]
		self.data = simulate(
			Seperation = False,
			Diffusion = False,
			Bloch = False,
			Pulse = False,
			ic = [str(_) for _ in self.ic],
			quiet = True)

	def test_rho1(self):
		self.assertEqual(self.data['rho1'][-1,:].tolist(), [self.ic[0]] * self.data['rho1'].shape[1])

	def test_rho2(self):
		self.assertEqual(self.data['rho2'][-1,:].tolist(), [self.ic[1]] * self.data['rho2'].shape[1])

	def test_rho3(self):
		self.assertEqual(self.data['rho3'][-1,:].tolist(), [self.ic[2]] * self.data['rho3'].shape[1])

class TestSSWithoutBloch(unittest.TestCase):

	@classmethod
	def setUpClass(self):
		self.ic = [0, 0.1, -0.2]
		self.data = simulate(
			Bloch = False,
			Pulse = False,
			ic = [str(_) for _ in self.ic],
			quiet = True)

	def test_rho1(self):
		self.assertEqual(self.data['rho1'][-1,:].tolist(), [self.ic[0]] * self.data['rho1'].shape[1])

	def test_rho2(self):
		self.assertEqual(self.data['rho2'][-1,:].tolist(), [self.ic[1]] * self.data['rho2'].shape[1])

	def test_rho3(self):
		self.assertEqual(self.data['rho3'][-1,:].tolist(), [self.ic[2]] * self.data['rho3'].shape[1])

class TestDiffusionWithoutBloch(unittest.TestCase):

	@classmethod
	def setUpClass(self):
		self.data = simulate(
			Bloch = False,
			Pulse = False,
			Seperation = False,
			ic = [
				'0',
				'-0.1',
				'0.1+0.1*cos(x[0]*3.14/7.512)'],
			quiet = True)

	def test_rho1(self):
		self.assertAlmostEqual(np.trapz(self.data['rho1'][0,:]), np.trapz(self.data['rho1'][-1,:]))

	def test_rho2(self):
		self.assertAlmostEqual(np.trapz(self.data['rho2'][0,:]), np.trapz(self.data['rho2'][-1,:]))

	def test_rho3(self):
		self.assertAlmostEqual(np.trapz(self.data['rho3'][0,:]), np.trapz(self.data['rho3'][-1,:]))

'''
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
'''

if __name__ == '__main__':
	unittest.main()
