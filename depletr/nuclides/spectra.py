"""
Spectra

Keep track of implemented nuclide libraries of various energy spectra
"""

from ._fast import get_fast_nuclib
from ._thermal import get_thermal_nuclib

#: dict of :code:`{"nuclib_name" : get_nuclib_function}`
SPECTRA = {"fast": get_fast_nuclib, "thermal": get_thermal_nuclib}
