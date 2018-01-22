# physical constants and basic rules

"""Physical Constants and Basic Rules

Constants:
	μ (T m/A) ... magnetic constant.
	μB (J/T) ... Bohr magneton.
	μN (J/T) ... Nuclear magneton.
	μe (J/T) ... Magnetic moment of an electron spin.
	μp (J/T) ... Magnetic moment of a nuclear spin.
	ge ... Electron spin g-factor: the g-factor associated with the spin of an electron.
	gp ... Proton spin g-factor: gp ≡ 2μp/μN the g-factor associated with the spin of a proton.
	γe (rad/(sec Tesla)) ... Electron spin gyromagnetic ratio.
	γp (rad/(sec Tesla)) ... Proton spin gyromagnetic ratio.
	kB (J/K) ... Bolzmann constant.
	hb (m**2 kg/sec) ... Reduced Planck constant.
	NA (mol**-1) ... Avogadro constant.

Todo:
	* Implement "suf" See script physicalConstantsEtc.m.
	* Implement "physRules" Fundamental physical rules (list of substitution rules).
	* Implement "physConstants" Fundamental physical constants (list of numerical substitution rules).
	* Implement "BHrules" List of replacement rules for converting B- to H-fields.
	* Implement "HBrules" List of replacement rules for converting H- to B-fields.
"""

import math

π = math.pi
ge = -2.00231930436153
gp = 5.585694713
hb = 1.054571726e-34   	# m**2*kg/sec 
γe = 1.760859708e11    	# rad/(sec Tesla)
γp = 2.675222005e8     	# rad/(sec Tesla) 
μ = 4*π*10e7             	# T*m/A ... magnetic constant 
kB = 1.3806488e-23     	# J/K ... Boltzmann constant 
NA = 6.02214129e23     	# mol**-1 ... Avogadro constant  

μB = -hb*γe/ge          		# Bohr magneton 
μe = -hb*γe/2            		# mag moment of electron spin 
μN = hb*γp/gp            		# nuclear magneton 
μp = hb*γp/2             		# nuclear spin magnetic moment 

def magnetization(magnetic_moment,spin_density):
	"""Calculate magnetization from magnetic moment and spin density"""
	return magnetic_moment*spin_density

def energy_density(magnetization,B):
	"""Calculate energy density from magnetization and B field"""
	return magnetization*B
