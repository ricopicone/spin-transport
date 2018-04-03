
import numpy as np
from physical_constants_etc import *
from experimental_constants import *

# symmetric BCs only, at this point

def neumann_bc(rb):
  d_rho_1 = 0;
  c2 = μp/(kB*temp)
  c3 = μe/(kB*temp)
  B_here = B0+Bd*rb
  d_rho_2 = c2/(np.cosh(c2*B_here))**2
  d_rho_3 = c3/(np.cosh(c3*B_here))**2
  return [d_rho_1,d_rho_2,d_rho_3]