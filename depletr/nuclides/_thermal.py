# Thermal Spectrum (PWR)

from .nuclide import Nuclide
from .nuclib import Nuclib
from .half_lives import *


def get_thermal_nuclib():
	"""Return the populated Nuclib for a thermal spectrum"""
	# Nuclide library to populate
	thermal = Nuclib("thermal", nu=2.5)
	# List to populate it with
	all_nuclides = []
	
	# Uranium
	u235 = Nuclide('U', 235)
	u235.sigma_f = 38.8
	u235.sigma_y = 8.7
	u235.nu = 2.4
	all_nuclides.append(u235)
	u238 = Nuclide('U', 238)
	u238.sigma_f = 0.103
	u238.sigma_y = 0.86
	all_nuclides.append(u238)
	u239 = Nuclide('u', 239)
	all_nuclides.append(u239)
	# Neptunium
	np237 = Nuclide('Np', 237)
	np237.sigma_f = 0.52
	np237.sigma_y = 33
	all_nuclides.append(np237)
	np238 = Nuclide('Np', 238)
	np238.sigma_f = 134
	np238.sigma_y = 13.6
	all_nuclides.append(np238)
	np239 = Nuclide('Np', 239)
	all_nuclides.append(np239)
	# Plutonium
	pu238 = Nuclide('Pu', 238)
	pu238.sigma_f = 2.4
	pu238.sigma_y = 27.7
	all_nuclides.append(pu238)
	pu239 = Nuclide('Pu', 239)
	pu239.sigma_f = 102.2
	pu239.sigma_y = 58.7
	pu239.nu = 2.9
	all_nuclides.append(pu239)
	pu240 = Nuclide('Pu', 240)
	pu240.sigma_f = 0.53
	pu240.sigma_y = 210.2
	all_nuclides.append(pu240)
	pu241 = Nuclide('Pu', 241)
	pu241.sigma_f = 102.2
	pu241.sigma_y = 40.9
	all_nuclides.append(pu241)
	pu242 = Nuclide('Pu', 242)
	pu242.sigma_f = 0.44
	pu242.sigma_y = 28.8
	all_nuclides.append(pu242)
	# Americium
	am241 = Nuclide('Am', 241)
	am241.sigma_f = 1.1
	am241.sigma_y = 110
	am241.metastable_branch_ratio = 0.11
	all_nuclides.append(am241)
	am242 = Nuclide('Am', 242)
	am242.sigma_f = 159
	am242.sigma_y = 301
	all_nuclides.append(am242)
	am242m = Nuclide('Am', 242)
	am242m.name += 'm'  # metastable
	am242m.latex += '$_m$'
	am242m.sigma_f = 595
	am242m.sigma_y = 137
	all_nuclides.append(am242m)
	am243 = Nuclide('Am', 243)
	am243.sigma_f = 0.44
	am243.sigma_y = 49
	all_nuclides.append(am243)
	# Curium
	cm242 = Nuclide('Cm', 242)
	cm242.sigma_f = 1.14
	cm242.sigma_y = 4.5
	all_nuclides.append(cm242)
	cm243 = Nuclide('Cm', 243)
	cm243.sigma_f = 88
	cm243.sigma_y = 14
	all_nuclides.append(cm243)
	cm244 = Nuclide('Cm', 244)
	cm244.sigma_f = 1.0
	cm244.sigma_y = 16
	all_nuclides.append(cm244)
	cm245 = Nuclide('Cm', 245)
	cm245.sigma_f = 116
	cm245.sigma_y = 17
	all_nuclides.append(cm245)
	
	for nuclide in all_nuclides:
		n = nuclide.name
		if n in ALPHA:
			nuclide.lambda_alpha = LN2/ALPHA[n]
		if n in BETAM:
			nuclide.lambda_betam = LN2/BETAM[n]
		if n in BETAP:
			nuclide.lambda_betap = LN2/BETAP[n]
		if n in GAMMA:
			nuclide.lambda_gamma = LN2/GAMMA[n]
		if nuclide.sigma_f and not nuclide.nu:
			nuclide.nu = thermal.nu
	
	# Finally, populate it.
	thermal.nuclides = all_nuclides
	return thermal
