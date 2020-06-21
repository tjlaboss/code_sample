"""
Nuclib

Module containing `Nuclib` class
"""

from collections import OrderedDict


class Nuclib:
	"""Nuclide library
	
	Add by doing :code:`Nuclib.nuclides = list_of_nuclides`.
	The `nuclides` dictionary will be built automatically.
	
	Parameters:
	-----------
	:param name:
		Brief name of the lib, such as "thermal"
	:type name: str
	
	:param nu:
		Average number of neutrons per fission
	:type nu: float
	
	Attributes:
	-----------
	:ivar nuclides:
		Dictionary of :code:`{"nuclide name" : Nuclide}`
	:vartype nuclides:
		collections.OrderedDict
	"""
	def __init__(self, name, nu):
		self.name = name
		self.nu = nu
		self._nuclides = None
	
	@property
	def nuclides(self):
		return self._nuclides
	
	@nuclides.setter
	def nuclides(self, list_of_nuclides):
		self._nuclides = OrderedDict()
		for nuc in list_of_nuclides:
			key = nuc.name.lower()
			self._nuclides[key] = nuc
	
	def __getattr__(self, item):
		try:
			return self._nuclides[item]
		except KeyError:
			errfmt = "Nuclib '{}' object has no attribute '{}'"
			raise AttributeError(errfmt.format(self.name, item))
		except TypeError:
			errfmt = "Nuclib '{}' has no nuclides yet."
			raise ValueError(errfmt.format(self.name))
