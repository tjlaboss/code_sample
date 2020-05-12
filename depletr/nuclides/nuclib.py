# Nuclib
#
# Nuclide libraries

from collections import OrderedDict


class Nuclib:
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
