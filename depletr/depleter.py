"""
Depleter

Module containing the `Depleter` class. It depletes nuclides.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import N_A
from scipy.linalg import expm as matrexp
from . import fuel, nuclides, plotter
from .dataset import DataSet


class Depleter:
	"""Depletion solver for the depletr module

	Parameters:
	-----------
	:type power: float
	:param power:
		Core power density in W/g
		
	:type enrichment: float
	:param enrichment:
		Initial U235 enrichment
	
	:type max_burnup: float
	:param max_burnup:
		Burnup limit to deplete to (i.e., EOC) in MW-d/kg-HM
	
	:type spectrum: str
	:param spectrum:
		"fast" or "thermal"
	
	:type mass: float
	:param mass:
		Heavy metal mass (kg) in the core at BOC
	
	Attributes:
	-----------
	:vartype fuel_type: str
	:ivar fuel_type:
		Hardcoded as "Uranium". Hope to implement other fresh fuels.
		
	:vartype data: nuclides.Nuclib
	:ivar data:
		Libary (either "fast" or "thermal") with nuclides defined
	"""
	def __init__(self, power, enrichment, max_burnup, spectrum, mass=1E8):
		self.power = power
		self.enrichment = enrichment / 100
		self.max_burnup = max_burnup
		self.mass = mass
		self.fuel_type = "Uranium"
		try:
			get_data = nuclides.SPECTRA[spectrum]
		except KeyError:
			errfmt = "Spectrum must be one of: {}"
			key_tuple = tuple(nuclides.SPECTRA.keys())
			raise ValueError(errfmt.format(key_tuple))
		self.data = get_data()
		self._scalar = mass*N_A/238*1E-24
		self._last_dataset = None
	
	
	def get_all_nuclides(self):
		"""Get all Nuclides in the system
		
		Returns:
		--------
		:rtype: dict_values of :class:`~nuclides.Nuclide`
		:return: All instances of Nuclide in the system
		"""
		return self.data.nuclides.values()
	
	
	def get_all_nuclide_names(self):
		"""Get the all nuclide names in the system
		
		Returns:
		--------
		:rtype: list of str
		:return: Names of all the Nuclide instances in the system
		"""
		return [n.name for n in self.get_all_nuclides()]
	
	
	def scale(self, quantities):
		"""Given the quantities of nuclides, find the actual concentrations
		
		One often goes back and forth between relative numbers (`quantities`)
		and absolute numbers (`concentrations`). This method serves to scale
		a vector of nuclide abundance--relative or absolute--to the number of
		nuclei in the reactor, assuming an average atomic mass of 238.
		
		Concentrations are scaled by 1E-24 (the ratio of 1 barn / 1 cm^2)
		for numerical precision reasons and easy micro/macro xs conversion.
		
		Parameter:
		----------
		:type quantities: np.ndarray of float; arbitrary units
		:param quantities:
			Vector of quantities of nuclides to scale to the reactor mass.
			Will be normalized: their units or absolute values do not matter.
		
		Returns:
		--------
		:rtype: np.ndarray of float; atoms * 1E-24
		:return:
			Vector of total number of nuclei
		"""
		return self._scalar*quantities/quantities.sum()
	
	
	def _deplete(self, nsteps, dt, ds, c0, fission_rate, verbose):
		"""Deplete the fuel at a constant fission rate
		
		Solves for the nuclide concentrations `nsteps` times, preserving
		a constant fission rate (equivalent to a constant power here).
		
		Parameters:
		-----------
		:type nsteps: int
		:param nsteps:
			Number of time steps to take per interval
		
		:type dt: float
		:param dt:
			Time step size, in seconds
		
		:type ds: depletr.DataSet
		:param ds:
			Data set to use for this depletion solution
		
		:type c0: np.ndarray, float, len=(ds.size)
		:param c0:
			Vector of absolute nuclide abundances at time zero.
			Scaled by factor of 1E-24.
		
		:type fission_rate: float
		:param fission_rate:
			Fission rate (fissions/second) to maintain at each solved time
		
		:type verbose: bool
		:param verbose
			If True, prints extra information
		
		Returns:
		--------
		concentrations: np.ndarray(dtype=float, shape=(ds.size, nsteps));
		                nuclei * 1E-24
			Nuclide concentrations at each time step.
		
		kinfvals: np.ndarray(dtype=float, len=nsteps)
			Vector of k-infinity at each time step.
		
		fluxvals: np.ndarray(dtype=float, len=nsteps); neutrons/cm^2/s
			Neutron flux at each time step.
		
		enrichvals: np.ndarray(dtype=float, len=nsteps)
			Enrichment (in atom fraction) of U235/(U235 + U238).
		"""
		m = ds.m
		l = ds.l
		kinfvals = np.zeros(nsteps)
		fluxvals = np.zeros(nsteps)
		enrichvals = np.zeros(nsteps)
		concentrations = np.zeros((ds.size, nsteps))
		enrichvals[0] = self.enrichment
		
		concentrations[:, 0] = c0
		fv = ds.get_fission_vector()
		nufv = ds.get_xs_vector("nu-fission")
		absv = ds.get_xs_vector("absorption")
		fission_xs = (c0*fv).sum()
		flux = fission_rate/fission_xs
		
		for k in range(nsteps):
			if k > 0:
				a = m*flux + l
				dn = matrexp(-a*dt)
				concentrations[:, k] = concentrations[:, k-1].dot(dn)
				enrichvals[k] = concentrations[0, k]/(concentrations[0, k] + concentrations[1, k])
			ck = concentrations[:, k]
			fission_xs = (ck*fv).sum()
			flux = fission_rate/fission_xs*1E-24
			fluxvals[k] = flux
			kinfvals[k] = (ck*nufv).sum()/(ck*absv).sum()
		
		if verbose:
			mass_reduct_238 = (concentrations[1, -1] - c0[1])/c0[1]
			print("U238 mass change: {:5.2%}".format(mass_reduct_238))
			mass_reduct_235 = (concentrations[0, -1] - c0[0])/c0[0]
			print("U235 mass change: {:5.2%}".format(mass_reduct_235))
			est235 = fuel.at_to_wt_uranium(enrichvals[-1])
			print("Final enrichment estimate: {:5.2%}".format(est235))
		return concentrations, kinfvals, fluxvals, enrichvals

	
	def deplete_fresh(self, nsteps, plots=1, verbose=True):
		"""Deplete fresh fuel uranium fuel to the burnup limit
		
		The nuclide list is initialized assuming pure U235/U238 fuel at
		Depleter.enrichment. Each step is depleted at a power of
		Depleter.power, which is used to calculate a constant fission rate
		and flux for each time step.
		
		Parameters:
		-----------
		:type nsteps: int
		:param nsteps:
			Number of time steps to use
		
		:type plots: {int | bool}, optional
		:param plots:
			Plotting mode:
			 -  0: (or anything Falsey) Do not make any plots
			 -  1: (or True) Make a single plot with 4 depletion subplots
			 -  2: For Debugging Only. Make matrix spy() plots and
			       initial heavy metal depletion plot.
		     - >2: Make 4 subplots
		    [Default: 1]
		
		:type verbose: bool (optional)
		:param verbose:
			Whether to print extra information during depletion.
			[Default: True]
		
		Returns:
		--------
		:rtype: np.ndarray, float; nuclides * 1E-24
		:returns:
			Nuclide concentrations at the end of the depletion (EOC).
		"""
		power_megawatt = self.power*self.mass*1E-6
		fission_rate = fuel.get_fission_rate(power_megawatt)
		all_nuclides = self.get_all_nuclides()
		quantities = fuel.give_me_fuel(self.fuel_type, self.enrichment, len(all_nuclides))
		time = fuel.give_me_fire(self.power, self.max_burnup)
		dt = time/nsteps
		
		ds = DataSet()
		ds.add_nuclides(all_nuclides, quantities)
		ds.build()
		c0 = self.scale(ds.get_initial_quantities())
		
		concs, kinf, flux, enrich = self._deplete(nsteps, dt, ds, c0, fission_rate, verbose)
		
		if plots:
			tvals = np.arange(0, nsteps*dt, dt)
			if plots == 1:
				fig = plt.figure()
				axa = fig.add_subplot(221)
				axf = fig.add_subplot(222)
				axk = fig.add_subplot(223)
				axu = fig.add_subplot(224)
			elif plots == 2:
				# Figure 1: Spy plots
				f = plt.figure()
				# Left: absorption matrix
				a = f.add_subplot(121)
				a.spy(ds.m.T)
				a.set_title("$[M]$\n")
				# Right: decay matrix
				b = f.add_subplot(122)
				b.spy(ds.l.T)
				b.set_title("$[L]$\n")
				# Figure 2: Actinide % of Initial Heavy Metal plot
				plt.figure()
				plotter.make_heavy_metal_plot(tvals, concs, all_nuclides, ax=None,
				                              deadend_actinides=True, fission_products=True,
				                              unit=(self.mass/self.power, "Burnup (MW-d/kg-HM)"))
				return concs[:, -1]
			else:
				axa = plt.figure().add_subplot(111)
				axf = plt.figure().add_subplot(111)
				axk = plt.figure().add_subplot(111)
				axu = plt.figure().add_subplot(111)
			# Top left: Actinide depletion
			plotter.make_actinides_plot(tvals, concs, all_nuclides, axa,
			                            fission_products=True, deadend_actinides=True)
			# Top Right: Enrichment and flux
			plotter.make_enrichment_flux_plot(tvals, enrich, flux, axf)
			# Bottom left: Approximate k-infinity
			plotter.make_kinf_plot(tvals, kinf, axk)
			# Bottom Right: Relative uranium depletion
			cvals = concs[0:2]
			nucs = (self.data.u235, self.data.u238)
			plotter.make_element_depletion_plot(tvals, cvals, nucs, axu,
			                                    relative=True, element="Uranium")
			if plots == 1:
				# Hide redundant xlabels to avoid crowding
				axa.set_xlabel("")
				axf.set_xlabel("")
		
		self._last_dataset = ds
		return concs[:, -1]
	
	
	def reprocess(self, quantities, which_elements, mox_frac):
		"""Reprocess spent into MOX
		
		Given the nuclide quantities in some spent fuel, take the ones we want
		and mix them with natural uranium to form mixed-oxide fuel (MOX).
		
		Parameters:
		-----------
		:type quantities: np.ndarray(dtype=float, len=NUM_NUCLIDES)
		:param quantities:
			Vector of nuclide quantities in the spent fuel,
			excluding lumped fission products and deadend actinides.
		
		:type which_elements: Iterable of str
		:param which_elements:
			Symbols for the elements to retain after chemical separation.
			For example, `which_elements=("Pu", "Np")` will retain all nuclides
			of the elements Plutonium and Neptunium.
		
		:type mox_frac: float
		:param mox_frac:
			Atom fraction of the retained elements to mix into the new fuel.
		
		Returns:
		--------
		:rtype: np.ndarray
		:return:
			Vector of nuclide abundances (*1E-24), including
			dead end actinides and lumped fission products.
		"""
		all_nuclides = self.get_all_nuclides()
		nuclides_names = self.get_all_nuclide_names()
		num = len(all_nuclides)
		mox = fuel.give_me_fuel("U", nuclides.constants.NATURAL_U235, num + 2)
		mox *= self._scalar
		indices = dict(zip(nuclides_names, range(num)))
		mox[indices["U235"]] *= (1 - mox_frac)
		mox[indices["U238"]] *= (1 - mox_frac)
		# New array with only the relevant nuclides.
		qnew = np.zeros_like(quantities)
		for nuc in all_nuclides:
			if nuc.element in which_elements:
				index = indices[nuc.name]
				qnew[index] = quantities[index]
		mox += mox_frac*self.scale(qnew)
		mox = self.scale(mox)
		return mox  # remember mox_frac is in atom fraction
	
	
	def _deplete_reloaded_fuel(self, quantities, nsteps, verbose):
		"""See `Depleter.reload()`."""
		if not self._last_dataset:
			self._last_dataset = DataSet()
			self._last_dataset.add_nuclides(self.get_all_nuclides(), quantities[:-2])
			self._last_dataset.build()
		ds = self._last_dataset
		power_megawatt = self.power*self.mass*1E-6
		fission_rate = fuel.get_fission_rate(power_megawatt)
		time = fuel.give_me_fire(self.power, self.max_burnup)
		dt = time/nsteps
		c0 = self.scale(quantities)
		return self._deplete(nsteps, dt, ds, c0, fission_rate, verbose)
	
	
	def reload(self, quantities, nsteps, verbose=False):
		"""Reload the core with known nuclides, then deplete to burnup limit
		
		When depleting a core that is not pure U235/U238, the nuclide
		abundance vector `quantities` must be provided. The ratio is
		preserved, while the magnitude is scaled to preserve the total
		heavy metal mass of fuel in the core (more or less).
		
		From then, the depletion calculation proceeds as in
		`Depleter.deplete_fresh()`.
		
		Parameters:
		-----------
		:type quantities: np.ndarray(dtype=float, len=NUM_NUCLIDES + 2)
		:param quantities:
			Vector of nuclides with which to fill the core.
			The last two entries are the "dead-end" (lumped) actinide
			and fission product, respectively.
			
		:type nsteps: int
		:param nsteps:
			Number of time steps to use
			
		:type verbose: bool (optional)
		:param verbose:
			Whether to print extra information during the depletion solution
			[Default: False]
		
		Returns;
		--------
		:rtype: np.ndarray(dtype=float, len=NUM_NUCLIDES + 2); nuclei * 1E-24
		:return:
			Vector of nuclide concentrations at the end of the depletion (EOC).
		"""
		concs, _, _, _ = self._deplete_reloaded_fuel(quantities, nsteps, verbose)
		return concs[:, -1]
	
	
	def reload_kinf(self, quantities, nsteps):
		"""`Depleter.reload()`, but returns k-infinities instead
		
		This function exists for optimization calculations.
		
		Returns;
		--------
		:rtype: np.ndarray(dtype=float, len=nsteps)
		:returns:
			Vector of k-infinity at every time step
		"""
		# For optimization purposes
		_, kinf, _, _ = self._deplete_reloaded_fuel(quantities, nsteps, verbose=False)
		return kinf
		
	
	def decay(self, quantities, nsteps, times):
		"""Let nuclides decay for a certain amount of time
		
		All nuclides will decay with their own half-lives,
		absent of any neutron flux.
		
		This method lazily uses an equal number of steps (`nsteps`) for every
		interval in `times`. It was sufficient for the 22.251 problem set,
		but should allow options to use an array of `nsteps` per interval,
		a fixed `dt`, or an array of `dt` per interval.
		
		Example:
		--------
		decayed = Depleter.decay(q_spent, 2, [10, 20, 100])
			Will deplete [0 -> 10], [10 -> 20], and [20 -> 100]
			with 2 steps: [5,  10,   15,   20,   60,  and  100]
		
		The resulting array `decayed` will represent nuclide concentrations
		at times 0, 10, 20, and 100.
		
		Parameters:
		-----------
		:type quantities: Iterable of float, len=(ALL_NUCLIDES + 2)
		:param quantities:
			Quantities of each Nuclide at the start of the decay period.
			The last two entries are the "dead-end" (lumped) actinide
			and fission product, respectively.
		
		:type nsteps: int
		:param nsteps:
			Number of time steps in each interval.
		
		:type times: np.ndarray of float, seconds
		:param times:
			Times to report Nuclide concentrations at.
			`nsteps` time steps will be taken between each time in the array.
		
		Returns:
		--------
		:rtype: np.ndarray of float, shape=(ALL_NUCLIDES + 1, len(times))
		:returns:
			Nuclide concentrations, including lumped actinides and fission products,
			at time 0 and each time requested in `times`.
		"""
		ds = DataSet()
		ds.add_nuclides(self.get_all_nuclides(), quantities[:-2])
		ds.build()
		# Ignore (lumped) fission products
		l = ds.l[:-1, :-1].copy()
		c0 = quantities[:-1].copy()
		# Deplete over the intervals of interest with nsteps each
		nt = len(times)
		nc = len(c0)
		results = np.zeros((nc, nt + 1))
		results[:, 0] = c0
		concentrations = np.zeros((nc, nt*nsteps + 1))
		concentrations[:, 0] = c0
		elapsed = 0
		i = 0
		for interval, time in enumerate(times):
			dt = (time - elapsed)/nsteps
			for k in range(nsteps):
				elapsed += dt
				dn = matrexp(-l*dt)
				concentrations[:, i+1] = concentrations[:, i].dot(dn)
			results[:, interval+1] = concentrations[:, i+1]
			i += 1
		return results
	
	
	def show(self):
		"""Shortcut for `matplotlib.pyplot.show()`  """
		return plt.show()
