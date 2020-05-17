# Fuel
#
# Functions related to fuel compositions


import numpy as np
from scipy.constants import eV
from .nuclides.half_lives import DAY


ELEMENTS = ("U", "Uranium")  # + ("Pu", "Plutonium")


def give_me_fuel(element, enrichment, num_nuclides):
	"""Get a nuclide vector for fresh fuel
	
	Parameters:
	-----------
	element: str
		Element name or atomic symbol.
		Only Uranium is currently supported.
	
	enrichment: float; wt%
		U235 enrichment
	
	num_nuclides: int
		Total number of nuclides to accommodate in the vector
	
	Returns:
	--------
	np.ndarray(dtype=float, len=num_nuclides)
		Array with nuclides populated for fresh fuel and 0 for everything else
	"""
	assert element in ELEMENTS, "Implemented fuel types: {}".format(ELEMENTS)
	assert num_nuclides >= 2, "You cannot have fewer than 2 nuclides."
	fuel_nuclides = np.zeros(num_nuclides)
	if element[0] == "U":
		at235 = wt_to_at_uranium(enrichment)
		at238 = 1 - at235
		fuel_nuclides[0] = at235
		fuel_nuclides[1] = at238
	else:
		# Plutonium, MOX, ...
		raise NotImplementedError(element)
	return fuel_nuclides


def give_me_fire(watts_per_gram, burnup):
	"""Get the time to achieve a required burnup
	
	Parameters:
	-----------
	watts_per_gram: float; W/g
		Power density of the reactor
	
	burnup: float; MW-d/kg-HM
		Desired fuel burnup at EOC
	
	Returns:
	--------
	float; seconds
		Time to reach `burnup` with a power density of `watts_per_gram`
	"""
	burnup *= DAY*1E3  # convert from MWd/kg to W-s/g
	return burnup/watts_per_gram


def give_me_that_which_i_desire():
	"""Ooh"""
	print("https://youtu.be/G1cjHbXdU0s?t=6")


def get_fission_rate(megawatts, mev_per_fission=200):
	"""Get the number of required fissions
	
	Parameters:
	-----------
	megawatts: float; MW
		Total reactor thermal power
	
	mev_per_fission: float; MeV/fission
		Average energy released per fission.
		TODO: Put this in the Nuclide class and data libraries.
		[Default: 200]
	
	Returns:
	--------
	float; MeV/s
		Fission rate required to achieve the desired power
	"""
	mevs = megawatts/eV
	return mevs/mev_per_fission


def wt_to_at_uranium(wt235):
	"""Convert weight fraction to atom fraction for uranium fuel
	
	Parameter:
	----------
	wt235: float
		Weight fraction U235
	
	Returns:
	--------
	float
		Atom fraction U235
	"""
	assert 0 <= wt235 <= 1, "Weight fraction must be on [0, 1]."
	wt238 = 1 - wt235
	# To be rigorous, replace with the actual atomic masses.
	at238 = wt238/238
	at235 = wt235/235
	return at235/(at235 + at238)


def at_to_wt_uranium(at235):
	"""Convert atom fraction to weight fraction for uranium fuel

	Parameter:
	----------
	at235: float
		Atom fraction U235

	Returns:
	--------
	float
		Weight fraction U235
	"""
	assert 0 <= at235 <= 1, "Atom fraction must be on [0, 1]."
	at238 = 1 - at235
	# To be rigorous, replace with the actual atomic masses.
	wt238 = at238*238
	wt235 = at235*235
	return wt235/(wt235 + wt238)
	

def wt_to_at_plutonium(wt239):
	return NotImplementedError("Plutonium fuel")


def at_to_wt_arbitrary(list_of_nuclides, list_of_ats):
	"""Convert multiple atom fractions to an array of weight fractions
	
	Parameters:
	-----------
	list_of_nuclides: Iterable of Nuclide
		Nuclides for which to convert atom fraction to weight fraction
	
	list_of_ats: Iterable of float
		Atom fractions of each nuclide
	
	Returns:
	--------
	np.ndarray(dtype=float)
		Array of weight fractions of each nuclide
	"""
	num = len(list_of_ats)
	total_at = sum(list_of_ats)
	list_of_wts = np.zeros(num)
	for i, (nuc, at) in enumerate(zip(list_of_nuclides, list_of_ats)):
		list_of_wts[i] = nuc.a*at/total_at
	list_of_wts /= list_of_wts.sum()
	return list_of_wts
