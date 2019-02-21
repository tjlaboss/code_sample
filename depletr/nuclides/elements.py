# Elements
#
# Convert between atomic number Z and element symbol

from collections import OrderedDict


SYMBOL = [None]*119
SYMBOL[0] = 'n'  # <neutron>
SYMBOL[1] = 'H'  # Hydrogen
SYMBOL[2] = 'He' # Helium
SYMBOL[3] = 'Li' # Lithium
SYMBOL[4] = 'Be' # Beryllium
SYMBOL[5] = 'B'  # Boron
SYMBOL[6] = 'C'  # Carbon
SYMBOL[7] = 'N'  # Nitrogen
SYMBOL[8] = 'O'  # Oxygen
SYMBOL[9] = 'F'  # Fluorine
SYMBOL[10] = 'Ne'  # Neon
SYMBOL[11] = 'Na'  # Sodium
SYMBOL[12] = 'Mg'  # Magnesium
SYMBOL[13] = 'Al'  # Aluminum
SYMBOL[14] = 'Si'  # Silicon
SYMBOL[15] = 'P'   # Phosphorus
SYMBOL[16] = 'S'   # Sulfur
SYMBOL[17] = 'Cl'  # Chlorine
SYMBOL[18] = 'Ar'  # Argon
SYMBOL[19] = 'K'   # Potassium
SYMBOL[20] = 'Ca'  # Calcium
SYMBOL[21] = 'Sc'  # Scandium
SYMBOL[22] = 'Ti'  # Titanium
SYMBOL[23] = 'V'   # Vanadium
SYMBOL[24] = 'Cr'  # Chromium
SYMBOL[25] = 'Mn'  # Manganese
SYMBOL[26] = 'Fe'  # Iron
SYMBOL[27] = 'Co'  # Cobalt
SYMBOL[28] = 'Ni'  # Nickel
SYMBOL[29] = 'Cu'  # Copper
SYMBOL[30] = 'Zn'  # Zinc
SYMBOL[31] = 'Ga'  # Gallium
SYMBOL[32] = 'Ge'  # Germanium
SYMBOL[33] = 'As'  # Arsenic
SYMBOL[34] = 'Se'  # Selenium
SYMBOL[35] = 'Br'  # Bromine
SYMBOL[36] = 'Kr'  # Krypton
SYMBOL[37] = 'Rb'  # Rubidium
SYMBOL[38] = 'Sr'  # Strontium
SYMBOL[39] = 'Y'   # Yttrium
SYMBOL[40] = 'Zr'  # Zirconium
SYMBOL[41] = 'Nb'  # Niobium
SYMBOL[42] = 'Mo'  # Molybdenum
SYMBOL[43] = 'Tc'  # Technetium
SYMBOL[44] = 'Ru'  # Ruthenium
SYMBOL[45] = 'Rh'  # Rhodium
SYMBOL[46] = 'Pd'  # Palladium
SYMBOL[47] = 'Ag'  # Silver
SYMBOL[48] = 'Cd'  # Cadmium
SYMBOL[49] = 'In'  # Indium
SYMBOL[50] = 'Sn'  # Tin
SYMBOL[51] = 'Sb'  # Antimony
SYMBOL[52] = 'Te'  # Tellurium
SYMBOL[53] = 'I'   # Iodine
SYMBOL[54] = 'Xe'  # Xenon
SYMBOL[55] = 'Cs'  # Cesium
SYMBOL[56] = 'Ba'  # Barium
SYMBOL[57] = 'La'  # Lanthanum
SYMBOL[58] = 'Ce'  # Cerium
SYMBOL[59] = 'Pr'  # Praseodymium
SYMBOL[60] = 'Nd'  # Neodymium
SYMBOL[61] = 'Pm'  # Promethium
SYMBOL[62] = 'Sm'  # Samarium
SYMBOL[63] = 'Eu'  # Europium
SYMBOL[64] = 'Gd'  # Gadolinium
SYMBOL[65] = 'Tb'  # Terbium
SYMBOL[66] = 'Dy'  # Dysprosium
SYMBOL[67] = 'Ho'  # Holmium
SYMBOL[68] = 'Er'  # Erbium
SYMBOL[69] = 'Tm'  # Thulium
SYMBOL[70] = 'Yb'  # Ytterbium
SYMBOL[71] = 'Lu'  # Lutetium
SYMBOL[72] = 'Hf'  # Hafnium
SYMBOL[73] = 'Ta'  # Tantalum
SYMBOL[74] = 'W'   # Tungsten
SYMBOL[75] = 'Re'  # Rhenium
SYMBOL[76] = 'Os'  # Osmium
SYMBOL[77] = 'Ir'  # Iridium
SYMBOL[78] = 'Pt'  # Platinum
SYMBOL[79] = 'Au'  # Gold
SYMBOL[80] = 'Hg'  # Mercury
SYMBOL[81] = 'Tl'  # Thallium
SYMBOL[82] = 'Pb'  # Lead
SYMBOL[83] = 'Bi'  # Bismuth
SYMBOL[84] = 'Po'  # Polonium
SYMBOL[85] = 'At'  # Astatine
SYMBOL[86] = 'Rn'  # Radon
SYMBOL[87] = 'Fr'  # Francium
SYMBOL[88] = 'Ra'  # Radium
SYMBOL[89] = 'Ac'  # Actinium
SYMBOL[90] = 'Th'  # Thorium
SYMBOL[91] = 'Pa'  # Protactinium
SYMBOL[92] = 'U'   # Uranium
SYMBOL[93] = 'Np'  # Neptunium
SYMBOL[94] = 'Pu'  # Plutonium
SYMBOL[95] = 'Am'  # Americium
SYMBOL[96] = 'Cm'  # Curium
SYMBOL[97] = 'Bk'  # Berkelium
SYMBOL[98] = 'Cf'  # Californium
SYMBOL[99] = 'Es'  # Einsteinium
SYMBOL[100] = 'Fm'  # Fermium
SYMBOL[101] = 'Md'  # Mendelevium
SYMBOL[102] = 'No'  # Nobelium
SYMBOL[103] = 'Lr'  # Lawrencium
SYMBOL[104] = 'Rf'  # Rutherfordium
SYMBOL[105] = 'Db'  # Dubnium
SYMBOL[106] = 'Sg'  # Seaborgium
SYMBOL[107] = 'Bh'  # Bohrium
SYMBOL[108] = 'Hs'  # Hassium
SYMBOL[109] = 'Mt'  # Meitnerium
SYMBOL[110] = 'Ds'  # Darmstadtium
SYMBOL[111] = 'Rg'  # Roentgenium
SYMBOL[112] = 'Cn'  # Copernicium
SYMBOL[113] = 'Nh'  # Nihonium
SYMBOL[114] = 'Fl'  # Flerovium
SYMBOL[115] = 'Mc'  # Moscovium
SYMBOL[116] = 'Lv'  # Livermorium
SYMBOL[117] = 'Ts'  # Go Vols
SYMBOL[118] = 'Og'  # Oganesson


Z = OrderedDict()
for zed, sym in enumerate(SYMBOL):
	Z[sym] = zed
