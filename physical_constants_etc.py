# physical constants and basic rules

'''
Old Mathematica usage messages ...
eventually, should turn these defs into classes or something

suf = "\nSee script physicalConstantsEtc.m."
physRules::usage = "Fundamental physical rules (list of substitution rules)." <> suf;
physConstants::usage = "Fundamental physical constants (list of numerical substitution rules)." <> suf;
BHrules::usage = "List of replacement rules for converting B- to H-fields" <> suf;
HBrules::usage = "List of replacement rules for converting H- to B-fields" <> suf;
μ::usage = "T m/A ... magnetic constant." <> suf;
μB::usage = "J/T ... Bohr magneton." <> suf;
μN::usage = "J/T ... Nuclear magneton." <> suf;
μe::usage = "J/T ... Magnetic moment of an electron spin." <> suf;
μp::usage = "J/T ... Magnetic moment of a nuclear spin." <> suf;
ge::usage = "Electron spin g-factor: the g-factor associated with the spin of an electron." <> suf;
gp::usage = "Proton spin g-factor: gp ≡ 2μp/μN the g-factor associated with the spin of a proton." <> suf;
γe::usage = "rad/(sec Tesla) ... Electron spin gyromagnetic ratio." <> suf;
γp::usage = "rad/(sec Tesla) ... Proton spin gyromagnetic ratio." <> suf;
kB::usage = "J/K ... Bolzmann constant." <> suf;
hb::usage = "m**2 kg/sec ... Reduced Planck constant." <> suf;
NA::usage = "mol**-1 ... Avogadro constant." <> suf;
'''

π = 3.14159265359
ge = -2.00231930436153
gp = 5.585694713
hb = 1.054571726 10**-34   	# m**2*kg/sec 
γe = 1.760859708 10**11    	# rad/(sec Tesla)
γp = 2.675222005 10**8     	# rad/(sec Tesla) 
μ = 4*π*10**-7             	# T*m/A ... magnetic constant 
kB = 1.3806488 10**-23     	# J/K ... Boltzmann constant 
NA = 6.02214129 10**23     	# mol**-1 ... Avogadro constant  

μB = -hb*γe/ge          		# Bohr magneton 
μe = -hb*γe/2            		# mag moment of electron spin 
μN = hb*γp/gp            		# nuclear magneton 
μp = hb*γp/2             		# nuclear spin magnetic moment 

def magnetization(magnetic_moment,spin_density):
	return magnetic_moment*spin_density

def energy_density(magnetization,B):
	return magnetization*B
