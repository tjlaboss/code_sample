"""
Nuclide

Module containing just the `Nuclide` class.
"""

import numpy as np
from . import elements
from . import constants


class Nuclide:
	"""A nuclide with several useful physical quantities

	Parameters:
	-----------
	:param element:
		atomic symbol for the element, e.g. "U"
	:type element: str
	
	:param a:
		mass number of the nucleus,    e.g. 235
	:type a: int

	Attributes:
	-----------
	:ivar name:           element+A, e.g., "U235" for Uranium-235
	:vartype name:        str
	:ivar alpha:          maximum fractional energy loss per collision
	:vartype alpha:       float
	:ivar xi:             mean logarithmic energy decrement
	:vartype xi:          float
	:ivar lambda:         decay constant (s^-1)
	:vartype lambda:      float
	:ivar halflife:       decay half-life (seconds)
	:vartype halflife:    float
	:ivar sigma_n:        (n, n) micro xs (barns)
	:vartype sigma_n:     float
	:ivar sigma_y:        (n, y) micro xs (barns)
	:vartype sigma_y:     float
	:ivar sigma_f:        (n, f) micro xs (barns)
	:vartype sigma_f:     float
	"""
	def __init__(self, element, a):
		element = element.title()
		self.element = element
		self.z = elements.Z[element]
		self.a = a
		self.name = element + str(a)
		self.latex = "${{}}^{{{}}}${}".format(a, element)
		# Scattering physics
		self._alpha = ((a - 1)/(a + 1))**2
		if a == 1:
			self._xi = 1
		else:
			self._xi = 1 + (self.alpha*np.log(self.alpha))/(1 - self.alpha)
		self.metastable_branch_ratio = 0
		# Cross sections
		self.sigma_n = 0  # scatter
		self.sigma_y = 0  # capture
		self.sigma_f = 0  # fission
		self.nu = 0       # neutrons/fission
		# Decay
		self.lambda_alpha = 0  # alpha decay
		self.lambda_betap = 0  # beta+ decay
		self.lambda_betam = 0  # beta- decay
		self.lambda_gamma = 0  # internal conversion
		self._lambda_total = None
		
	def __str__(self):
		return "Nuclide: {} (kinf = {:5.2f})".format(self.name, self.kinf)
	
	@property
	def alpha(self):
		return self._alpha
	
	@property
	def xi(self):
		return self._xi
	
	@property
	def nu_sigma_f(self):
		return self.nu*self.sigma_f
	
	@property
	def sigma_a(self):
		return self.sigma_f + self.sigma_y
	
	@property
	def kinf(self):
		return self.nu_sigma_f / self.sigma_a if self.sigma_a else 0
	
	@property
	def lambda_total(self):
		if self._lambda_total is None:
			return self.lambda_betam + self.lambda_betap + \
			       self.lambda_alpha + self.lambda_gamma
	
	@lambda_total.setter
	def lambda_total(self, l):
		self._lambda_total = l
	
	@property
	def halflife(self):
		return np.log(2)/self._lambda_total
	
	@halflife.setter
	def halflife(self, t12):
		self._lambda_total = np.log(2)/t12

	def capture(self):
		"""Get the daughter nuclides from a neutron capture
		
		For nuclei where we model excited but metastable states, there exists
		a branch ratio for producing the ground state or the metastable state.
		
		Returns:
		--------
		:return: Names and branch ratios of the daughter nuclide at the ground
		         and metastable (excited) states, such as:
		         ((name_g, branch_ratio_g), (name_m, branch_ratio_m))
		:rtype:  Nested tuples of the structure: ((str, str), (str, str))
		"""
		daughter = self.element + str(self.a + 1)
		state0 = (daughter, 1 - self.metastable_branch_ratio)
		statem = (daughter + 'm', self.metastable_branch_ratio)
		return state0, statem
	
	def decay_betap(self):
		"""Get the daughter nuclide from a Beta+ decay
		
		Returns:
		--------
		:return: Name of the daughter nuclide
		:rtype:  str or Nonetype
		"""
		ep = elements.SYMBOL[self.z - 1]
		return ep + str(self.a)
	
	def decay_betam(self):
		"""Get the daughter nuclide from a Beta- decay
		
		Returns:
		--------
		:return: Name of the daughter nuclide
		:rtype:  str or Nonetype
		"""
		em = elements.SYMBOL[self.z + 1]
		return em + str(self.a)
	
	def decay_alpha(self):
		"""Get the daughter nuclide from an alpha decay
		
		Returns:
		--------
		:return: Name of the daughter nuclide
		:rtype:  str or Nonetype
		"""
		e = elements.SYMBOL[self.z - 2]
		a = self.a - 4
		return e + str(a)
	
	def decay_gamma(self):
		"""Get the daughter nuclide from internal conversion or gamma decay
		
		Returns:
		--------
		:return: Name of the daughter nuclide
		:rtype:  str or Nonetype
		"""
		if self.name[-1] == "m":
			return self.name[:-1]
	
	def fission(self):
		"""Get the fission product, if fissionable.
		
		This returns a list instead of a single string, in case we want
		to return 2 (or 3) fission products.
		
		Returns:
		--------
		:return: List containing just the lumped fission product
		:rtype:  list of str, or Nonetype
		"""
		if self.sigma_f:
			return [constants.FISSION_PRODUCT]
	
	def get_all_daughters(self):
		"""Nuclides, lock up your daughters
		
		Returns:
		--------
		:return:  daughter nuclides from all possible decays
		:rtype:   list of str
		"""
		daughters = []
		if self.lambda_alpha:
			daughters.append(self.decay_alpha())
		if self.lambda_betam:
			daughters.append(self.decay_betam())
		if self.lambda_betap:
			daughters.append(self.decay_betap())
		if self.lambda_gamma:
			daughters.append(self.decay_gamma())
		return daughters
