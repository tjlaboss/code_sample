# Data Set

from collections import OrderedDict
import numpy as np
from . import nuclides


class DataSet:
	"""Depletion data collection
	
	Populate the DataSet using the method `DataSet.add_nuclide()` or
	`DataSet.add_nuclides()`. Once all nuclides are in, perform the
	necessary nuclide linking using `DataSet.build()`.
	
	These attributes are only available after the DataSet has been built.
	
	Attributes:
	-----------
	size: int
		Size of any axis of the DataSet's 'M' and 'L' arrays.
		This is equal to the number of nuclides in the system plus two
		(one for lumped actinides and one for lumped fission products).
	
	m: np.ndarray, float (barns)
		'M' matrix of microscopic cross sections. Shape: (size, size)
	
	l: np.ndarray, float (s^-1)
		'L' matrix of decay constants (lambdas). Shape: (size, size)
	
	"""
	def __init__(self):
		self._nuclides = OrderedDict()
		self._q = OrderedDict()  # quantities
		self._built = False
		# Variables set after building
		self._size = None
		self._m = None
		self._l = None
		self._rxn_vectors = {}
	
	def _check_built(self):
		assert self._built, "You must build this DataSet first."
	
	@property
	def size(self):
		self._check_built()
		return self._size
	
	@property
	def m(self):
		self._check_built()
		return self._m
	
	@property
	def l(self):
		self._check_built()
		return self._l
	
	def add_nuclide(self, nuc, q=0):
		"""Add a nuclide to the set
		
		Parameters:
		-----------
		nuc: Nuclide
			Nuclide to add to the system. If it is already present,
			an error will be raised.
			
		q: float; optional
			Quantity of the nuclide to add. Arbitrary units.
			Specify a quantity of 0 to add the nuclide to the depletion chain.
			[Default: 0]
		
		"""
		if self._built:
			raise ValueError("Cannot add nuclides; matrix already built.")
		key = nuc.name
		if key in self._nuclides:
			errstr = "Nuclide {} already exists.".format(key)
			raise ValueError(errstr)
		self._nuclides[key] = nuc
		self._q[key] = q
	
	def add_nuclides(self, list_of_nuclides, list_of_quantities=None):
		"""Add multiple nuclides to the set
		
		Parameters:
		-----------
		list_of_nuclides: Iterable of Nuclide
			Nuclides to add to the system. If any are already present,
			an error will be raised.
		list_of_quantities: Iterable of float; optional
			Quantities of each nuclide in `list_of_nuclides.`
			Specify quantities of zero to add nuclides to the depletion chain.
			[Default: None --> all Nuclides at 0]
		
		"""
		n_nuc = len(list_of_nuclides)
		if list_of_quantities is None:
			list_of_quantities = [0]*n_nuc
		n_qnt = len(list_of_quantities)
		errstr = "Cannot add_nuclides(); nuclide list has {n} entries while quantity list has {q} entries."
		assert n_nuc == n_qnt, errstr.format(n=n_nuc, q=n_qnt)
		for nuc, q in zip(list_of_nuclides, list_of_quantities):
			self.add_nuclide(nuc, q)
		
	def build(self):
		"""Build the depletion matrices M and L."""
		indices = dict(zip(self._nuclides.keys(), range(len(self._nuclides))))
		indices[nuclides.FISSION_PRODUCT] = -1
		indices[nuclides.DEADEND_ACTINIDE] = -2
		
		n = len(indices)  # including the lumped boys
		M = np.zeros((n, n))  # sigma
		L = np.zeros((n, n))  # lambda
		
		warnstr = "{} daughter {} of nuclide {} not in data set."
		for i, (key, nuclide) in enumerate(self._nuclides.items()):
			# Populate the diagonal with the absorption xs
			M[i, i] = nuclide.sigma_a
			# Include the capture daughters...
			for daughter, branch_ratio in nuclide.capture():
				if not branch_ratio:
					continue
				if daughter in self._nuclides:
					j = indices[daughter]
					M[i, j] = -nuclide.sigma_y*branch_ratio
				elif nuclide.sigma_y:
					print(warnstr.format("capture", daughter, nuclide.name))
					j = indices[nuclides.DEADEND_ACTINIDE]
					M[i, j] = -nuclide.sigma_y*branch_ratio
			# ...and the fission products.
			if nuclide.sigma_f:
				j = indices[nuclides.FISSION_PRODUCT]
				M[i, j] = -nuclide.sigma_f
			
			# Populate the decay matrix with the lambdas
			L[i, i] = nuclide.lambda_total
			if nuclide.lambda_betam:
				daughter = nuclide.decay_betam()
				if daughter in self._nuclides:
					j = indices[daughter]
				else:
					print(warnstr.format("Beta-", daughter, nuclide.name))
					j = indices[nuclides.DEADEND_ACTINIDE]
				L[i, j] += -nuclide.lambda_betam
			if nuclide.lambda_betap:
				daughter = nuclide.decay_betap()
				if daughter in self._nuclides:
					j = indices[daughter]
				else:
					print(warnstr.format("Beta+", daughter, nuclide.name))
					j = indices[nuclides.DEADEND_ACTINIDE]
				L[i, j] += -nuclide.lambda_betap
			if nuclide.lambda_alpha:
				daughter = nuclide.decay_alpha()
				if daughter in self._nuclides:
					j = indices[daughter]
				else:
					print(warnstr.format("alpha", daughter, nuclide.name))
					j = indices[nuclides.DEADEND_ACTINIDE]
				L[i, j] += -nuclide.lambda_alpha
			if nuclide.lambda_gamma:
				# e.g. if it's Am-242m
				daughter = nuclide.decay_gamma()
				if daughter in self._nuclides:
					j = indices[daughter]
				else:
					j = indices[nuclides.DEADEND_ACTINIDE]
				L[i, j] += -nuclide.lambda_gamma
		
		self._built = True
		self._size = n
		self._m = M
		self._l = L
	
	def get_initial_quantities(self):
		"""Get the vector of quantities at time zero.
		
		Returns:
		--------
		None, or np.ndarray of float
			Returns None if the DataSet has not been built.
			Otherwise, returns a 1D array of the initial quantities.
		"""
		if self._built:
			if nuclides.DEADEND_ACTINIDE not in self._nuclides:
				self._q[nuclides.DEADEND_ACTINIDE] = 0
			if nuclides.FISSION_PRODUCT not in self._nuclides:
				self._q[nuclides.FISSION_PRODUCT] = 0
			vals = tuple(self._q.values())
			return np.array(vals)
	
	def _build_xs_vector(self, rxn):
		"""Build the vector of all nuclides' cross sections for a reaction
		
		Parameter:
		----------
		rxn: str
			Reaction to find. One of: {"fission" | "nu-fission" | "absorption" | "capture"}
		
		Returns:
		--------
		np.ndarray of float
			Array of cross sections of `rxn` for each nuclide.
		"""
		vector = np.zeros(self._size)
		for i, nuclide in enumerate(self._nuclides.values()):
			if rxn == "fission":
				xs = nuclide.sigma_f
			elif rxn == "nu-fission":
				xs = nuclide.nu_sigma_f
			elif rxn == "absorption":
				xs = nuclide.sigma_a
			elif rxn == "capture":
				xs = nuclide.sigma_y
			else:
				raise NotImplementedError(rxn)
			vector[i] = xs
		# Leave deadend actinides and fission products at 0.
		return vector
	
	def get_xs_vector(self, rxn):
		"""Get the vector of all nuclides' cross sections for a reaction
		
		If the vector exists, it will be looked up.
		Otherwise, it will be created.
		
		Parameter:
		----------
		rxn: str
			Reaction to find. One of: {"fission" | "nu-fission" | "absorption" | "capture"}
		
		Returns:
		--------
		np.ndarray of float
			Array of cross sections of `rxn` for each nuclide.
		"""
		self._check_built()
		if rxn not in self._rxn_vectors:
			self._rxn_vectors[rxn] = self._build_xs_vector(rxn)
		return self._rxn_vectors[rxn]
	
	def get_fission_vector(self):
		"""Get the vector of all nuclides' fission cross sections"""
		return self.get_xs_vector("fission")
