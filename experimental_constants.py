# define experimental constants
#   UW experiment parameters

from physical_constants_etc import *

"""
below are the legacy Mathematica usage messages
suf = "\nSee script: experimentalConstants.m"
experimentalRules::usage = "List of replacement rules for experimental parameters." <> suf;
B0::usage = "tesla ... Background field Dougherty2000." <> suf;
grad::usage = "tesla/meter ... Field gradient Dougherty2000." <> suf;
g1::usage = "dimensionless ... Nuclear spin g-factor." <> suf;
Γ1::usage = "meter**2/sec ... Dipole energy transport coefficient." <> suf;
Γ2::usage = "meter**2/sec ... Nuclear polarization transport coefficient." <> suf;
Γ3::usage = "meter**2/sec ... Electron polarization transport coefficient." <> suf;
ΓOZ::usage = "meter**2/sec ... Single transport coefficient of the OZ-Ansatz. Used when approximately true (e.g. single spin-species)." <> suf;
Γp::usage = "meter**2/sec ... Nuclear polarization diffusion coefficient (Dougherty2000)." <> suf;
Γe::usage = "meter**2/sec ... Electron polarization diffusion coefficient (Dougherty2000)." <> suf;
Bd2::usage = "tesla ... Mean dipole field from nuclear spins." <> suf;
Bd3::usage = "tesla ... Mean dipole field from electron spins." <> suf;
Bd::usage = "tesla ... Mean total dipole field from all spin species." <> suf;
temp::usage = "kelvin ... Experimental (physical) temperature." <> suf;
MwPS::usage = "g/mol ... Molar mass of polystyrene (not using polymer chain)." <> suf;
dPS::usage = "g/mL ... Volumetric density of polystyrene (from bottle)" <> suf;
nAMPS::usage = "Number of **1H atoms per molecule of polystyrene (not considering polymer chain)." <> suf;
concPS::usage = "Concentration of polystyrene in solution." <> suf;
δ2::usage = "1/m**3 ... Nuclear spin density." <> suf;
MwDPPH::usage = "g/mol ... Molar mass of DPPH (Wikipedia)." <> suf;
dPS::usage = "g/mL (g/cm**3) ... Volumetric density of DPPH (Wikipedia)" <> suf;
nAMDPPH::usage = "Number of free radicals per molecule of DPPH." <> suf;
concDPPH::usage = "Concentration of DPPH in polystyrene solution." <> suf;
δ3::usage = "1/m**3 ... Electron spin density." <> suf;
"""

# general
temp = 10 # Kelvin ... physical temperature of sample 

# electron spin concentration
MwDPPH = 394.32 # g/mol ... molar mass of DPPH (wikipedia) 
dDPPH = 1.4 # g/cm**3 ... density of DPPH (wikipedia) 
nAMDPPH = 1 # just one free radical per molecule 
concDPPH = .04 # concentration DPPH 
den3 = 1/MwDPPH*dDPPH*1e6*NA*nAMDPPH # 1/m**3 
δ3 = concDPPH*den3 # 1/m**3 

# nuclear spin concentration
#MwPS = 14 # amu ... molecular mass of polystyrene 
MwPS = ( # g/mol ... molar mass of polystyrene assuming C8H8 (Wikipedia)
    12*0 +  # no spin-1/2 from C atoms  + 
    1*8)    # H atoms  
dPS = 1.047 # g/mL ... density of polystyrene (from bottle) 
nAMPS = 2 # number of **1H atoms per molecule for polystyrene 
concPS = 1 - concDPPH           # concentration polystyrene 
den2 = 1/MwPS*dPS*1e6*NA*nAMPS # 1/m**3 
δ2 = concPS*den2 # 1/m**3 
# δ2 = 1.34*10**26 1/m**3    # Dougherty2000  
δb = δ3/δ2 # dimensionless spin concentration ratio

# identifying spins
γ2 = γp # species 2 is proton
γ3 = γe # species 3 is electron
γb = γ3/γ2 # dimensionless gyromagnetic ratio
g1 = 2.79 # dimensionless nuclear spin g-factor 

# magnetic fields
B0 = 0.0893 # T ... background field
grad = 44e3 # T/m ... background field gradient
Bd2 = μ/(4 * π)*hb*γ2*δ2 # T ... nuclear dipole field
Bd3 = μ/(4 * π)*hb*γ3*δ3 # T ... electron dipole field Dougherty2000 
Bd = Bd2 + Bd3 # T ... not the worst guess
Bb = 1 # dimensionless B ratio ... grad/Bd /(grad/Bd)

# transport rates
Γp = μ/(4 * π)*hb*γ2**2*δ2**(1/3)
Γe = μ/(4 * π)*hb*γ3**2*δ3**(1/3)
Γ2 = Γp
Γ3 = Γe
Γb = Γ3/Γ2 # dimensionless transport coefficient ratio
ΓOZ = Γ2 # used in case of single transport rate

# c ratio
cb = Bb*(1+δb)/(1+γb*δb)

# relaxation rates
T1e = 30.3e-6 # sec ... electron T1 at 10K (Dougherty2000)
T2e = 20e-9 # sec ... electron T2 at 10K (Dougherty2000)
T1p = .1 # sec ... nuclear T1
T2p = T2e # sec ... nuclear T2 ... TODO look up good value
T12 = T1p
T22 = T2p
T13 = T1e
T23 = T2e

def nondimensionalize_time(time_variable):
    """returns dimensionless time"""
    return Γ2*(grad/Bd)**2*time_variable

def nondimensionalize_space(space_variable):
    """returns dimensionless space"""
    return grad/Bd*space_variable
