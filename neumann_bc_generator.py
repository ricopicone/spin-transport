
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

if __name__ == '__main__':
  rb = float(input('Enter r bar: '))
  print(neumann_bc(rb))
