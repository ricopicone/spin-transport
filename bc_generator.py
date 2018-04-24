import sys
import numpy as np
from physical_constants_etc import *
from experimental_constants import *

# symmetric BCs only, at this point

def neumann_bc(rb):
  d_rho_1 = 0;
  c2 = μp/(kB*temp)
  c3 = μe/(kB*temp)
  B_here = B0+Bd*rb
  d_rho_2 = Bd*c2/(np.cosh(c2*B_here))**2
  d_rho_3 = Bd*c3/(np.cosh(c3*B_here))**2
  #print("%e %e %e" % (d_rho_1,d_rho_2,d_rho_3))
  return [d_rho_1,d_rho_2,d_rho_3]

def dirichlet_bc(rb):
  c2 = μp/(kB*temp)
  c3 = μe/(kB*temp)
  B_here = B0+Bd*rb
  # Langevin constants
  Ω_10 = (1+γb*δb)*(μp*Bd)/(kB*temp)
  Ω_20 = -μp*B0/(kB*temp)
  Ω_30 = -δb*μe*B0/(kB*temp)
  # corresponding rhos
  rho_1 = -np.tanh(Ω_10/(1+δb))
  rho_2 = -np.tanh(Ω_20 - 1/(1+γb*δb) * (B_here - B0)/Bd * Ω_10)
  rho_3 = -np.tanh(Ω_30/δb - γb*δb/(1+γb*δb) * (B_here - B0)/Bd * Ω_10)
  #print("%e %e %e" % (d_rho_1,d_rho_2,d_rho_3))
  return [rho_1,rho_2,rho_3]

if __name__ == '__main__':
  rb = float(input('Enter r bar: '))
  try:
    if sys.argv[1].lower() == 'neumann':
      print(neumann_bc(rb))
    if sys.argv[1].lower() == 'dirichlet':
      print(dirichlet_bc(rb))
    else:
      print('{} is not an allowed boundary condition type, Use Dirichlet or Neumann'.format(sys.argv[1]))
  except IndexError:
    print('You must specicify either Neumann or Dirichlet BC as the first command line argument')
