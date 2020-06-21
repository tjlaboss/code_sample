"""
Printer

Printing functions for depletr
"""

import numpy as np


def print_decay_results(names, times, concentrations, other_actinides="other", **kwargs):
	"""Print a fancy table of nuclide inventory decay
	
	Parameters:
	-----------
	:type names: Iterable of str
	:param names:
		Names of each nuclide type in the repository
	
	:type times: Iterable of float; years
	:param times:
		Years after shutdown for each column.
	
	:type concentrations:
		np.ndarray(dtype=float, shape=[len(names)+1, len(times)+1]
	:param concentrations:
		Nuclide concentrations over time.
		Includes dead-end actinides. Excludes lumped fission products.
	
	:type other_actinides: str; optional
	:param other_actinides:
		What to call the lumped dead-end actinide in the table. Keep it short.
		[Default: "other"]
	
	:type kwargs: dict
	:param kwargs:
		Keyword arguments to pass down to the print() function.
	"""
	assert concentrations.shape == (len(names)+1, len(times)+1), \
		"Misshapen data. Format is: concentrations[nuclides+1, times+1]"
	names += [other_actinides]
	ndec = 16
	nchar = len(max(names, key=len))
	print("\n\nDecay after __ years:\n", **kwargs)
	
	fmt_header = '{:^' + str(ndec) + '}|'
	headstr = " "*nchar + " | "
	for t in [0] + list(times):
		headstr += fmt_header.format(t)
	print(headstr, **kwargs)
	print("-"*len(headstr), **kwargs)
	datastr = np.array2string(concentrations, max_line_width=np.inf, separator=' | ')
	datastr = datastr.replace('[', ' ').replace(']', ' ').split('\n')
	fmt_nuclide = '{:^' + str(nchar) + '} |'
	for i, name in enumerate(names):
		print(fmt_nuclide.format(name) + datastr[i], **kwargs)
	print(**kwargs)
