# absolute module imports
from mbfit.exceptions import InvalidValueError

"""
Basic unit constants
"""
avogadro = 6.02214129E+23    # NIST 2010
au_to_ev = 27.21138505       # NIST 2010
au_to_joule = 4.35974434E-18 # NIST 2010
cal_to_joule = 4.184
bohr = 0.52917721067e-10 # meter, 2014 CODATA
bohr_to_ang = bohr * 1.0e10 # 0.52 ang / 1 bohr
ang_to_bohr = 1/bohr_to_ang # 1.82 bohr / 1 ang
kelvin_to_au =  3.16683 * 1e-6 

autocm = 2.194746313e5 # Hartree to wavenumbers (cm-1) - TODO: derive from fundamental constants
cmtoau = 4.5563352527e-6 # wavenumbers (cm-1) to Hartree
mass_electron_per_mass_proton = 1822.88839 # mass of proton / mass of electron
"""
Constants derived from basic unit constants
"""                                                                      
au_to_kcal = au_to_joule * avogadro * 1E-03 / cal_to_joule
ev_to_kcal = au_to_kcal / au_to_ev
au_times_bohr6_to_kcal_times_ang6 = au_to_kcal * (bohr_to_ang ** 6) # Converts au*(bohr^6) to (kcal/mol)*(angstrom^6)

lone_pairs = ["X", "Y", "Z"]

# list of atomic symbols listed in order of atomic number, 0th item has atomic number 1, 1st item has atomic number 2,
# etc.
atomic_symbols = [
    "H",                                                                                                                                    "He",
    "Li",   "Be",                                                                                   "B",    "C",    "N",    "O",    "F",    "Ne",
    "Na",   "Mg",                                                                                   "Al",   "Si",   "P",    "S",    "Cl",   "Ar",
    "K",    "Ca",   "Sc",   "Ti",   "V",    "Cr",   "Mn",   "Fe",   "Co",   "Ni",   "Cu",   "Zn",   "Ga",   "Ge",   "As",   "Se",   "Br",   "Kr",
    "Rb",   "Sr",   "Y",    "Zr",   "Nb",   "Mo",   "Tc",   "Ru",   "Rh",   "Pd",   "Ag",   "Cd",   "In",   "Sn",   "Sb",   "Te",   "I",    "Xe",
    "Cs"
]

# list of atomic masses listed in order of atomic number
atomic_masses = [ # Source ptable.com (not a final source)
    1.008,                                                                                                                                  4.0026,
    6.94,   9.0122,                                                                                 10.81,  12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305,                                                                                 26.982, 28.085, 30.974, 32.06,  35.45,  39.948,
    39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38,  69.723, 72.630, 74.922, 78.971, 79.904, 83.798,
    85.468, 87.62,  88.906, 91.224, 92.906, 95.95,  98,     101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60, 126.90, 131.29,
    132.91
]

# list of atomic radii (angstroms) in order of atomic number
atomic_radii = [ # Source ptable.com
    0.53,                                                                                                                                   0.31,
    1.67,   1.12,                                                                                   0.87,   0.67,   0.56,   0.48,   0.42,   0.38,
    1.90,   1.45,                                                                                   1.18,   1.11,   0.98,   0.88,   0.79,   0.71,
    2.43,   1.94,   1.84,   1.76,   1.71,   1.66,   1.61,   1.56,   1.52,   1.49,   1.45,   1.42,   1.36,   1.25,   1.14,   1.03,   0.94,   0.88,
    2.65,   2.19,   2.12,   2.06,   1.98,   1.90,   1.93,   1.78,   1.73,   1.69,   1.65,   1.61,   1.56,   1.45,   1.33,   1.23,   1.15,   1.08,
    2.98
]

# list of covalent radii (angstroms) in order of atomic number
covalent_radii = [ # Source ptable.com
    0.37,                                                                                                                                   0.32,
    1.34,   0.90,                                                                                   0.82,   0.77,   0.75,   0.73,   0.71,   0.69,
    1.54,   1.30,                                                                                   1.18,   1.11,   1.06,   1.02,   0.99,   0.97,
    1.96,   1.74,   1.44,   1.36,   1.25,   1.27,   1.39,   1.25,   1.26,   1.21,   1.38,   1.31,   1.26,   1.22,   1.19,   1.16,   1.14,   1.10,
    2.11,   1.92,   1.62,   1.48,   1.37,   1.45,   1.56,   1.26,   1.35,   1.31,   1.53,   1.48,   1.44,   1.41,   1.38,   1.35,   1.33,   1.30,
    2.25   
]

# list of van der Waals radii (vdw) (angstroms) in order of atomic number
vdw_radii = [ # Source: ptable.com for elements unless in list below.
			  # Be, B, Al, Ca (Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/)
        # Rb, Cs (Source: doi:10.1021/jp8111556)
			  # vdw_radius is -1.0, if a definite value cannot be found. 
    1.2,                                                                                                                                    1.4,
    1.82,   1.53,                                                                                   1.92,   1.70,   1.55,   1.52,   1.47,   1.54,
    2.27,   1.73,                                                                                   1.84,   2.10,   1.80,   1.80,   1.75,   1.88,
    2.75,   2.31,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   1.63,   1.40,   1.39,   1.87,   -1.0,   1.85,   1.90,   1.85,   2.02,
    3.03,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   1.63,   1.72,   1.58,   1.93,   2.17,   -1.0,   2.06,   1.98,   2.16,
    3.43
]

# list of free atomic polarizabilities (bohr^3) in order of atomic number
free_polarizabilities = [ # Source: 2018 Table of static dipole polarizabilities of the neutral elements in the periodic table (Peter Schwerdtfeger & Jeffrey K. Nagle)
                          # DOI: 10.1080/00268976.2018.1535143
                          # https://www.tandfonline.com/doi/pdf/10.1080/00268976.2018.1535143?needAccess=true
     4.50711,                                                                                                              1.38375,
    164.1125, 37.74,                                                                        20.5, 11.3,  7.4,   5.3, 3.74, 2.66110,
       162.7,  71.2,                                                                        57.8, 37.3,   25,  19.4, 14.6,  11.083,
       289.8, 160.8,    97,   100,    87,    83,    68,    62,    55,    49,  46.5, 38.67,    50,   40,   30,  28.9,   21,   16.78,
       319.8, 197.2,   162,   112,    98,    87,    79,    72,    66, 26.14,    55,    46,    65,   53,   43,    38, 32.9,   27.32,
       400.9
]

# list of free atomic polarizabilities (angstrom^3) in order of atomic number
# Source: calculated using ccsd(t) with MolPRO (might not be accurately sourced, defer to Marc Riera)
ccsdt_free_polarizabilities = {
    "H":    0.66582,
    "B":    3.64084,
    "C":    1.84613,
    "N":    1.08053,
    "O":    0.86381,
    "F":    0.50682,
    "P":    3.72507,
    "S":    3.24768
}

def symbol_to_number(symbol):
    """
    Converts an atomic symbol to an atomic number.

    Args:
        symbol              - The 1 or 2 letter atomic symbol to convert to an atomic number. For example: "He", "F".
                Case non-sensitive.

    Returns:
        The atomic number for the atom specified by the given symbol.
    """

    # Make first letter uppercase and second letter (if exists) lowercase
    symbol = symbol[:1].upper() + symbol[1:].lower()

    try:
        # convert atomic symbol to number by performing a look-up in the list of atomic symbols
        return atomic_symbols.index(symbol) + 1
    except ValueError:
        # if the given symbol was not found in the list, then it is an invalid symbol
        raise InvalidValueError("atomic symbol", symbol, "a valid 1 or 2 letter atomic symbol") from None

def number_to_symbol(number):
    """
    Converts an atomic number to an atomic symbol.

    Args:
        number              - The atomic number to convert to an atomic symbol.

    Returns:
        The atomic symbol of the atom with the given atomic number.
    """

    # if the number is 0 or less, it is an invalid atomic number
    if number < 1:
        raise InvalidValueError("atomic number", number, "greater than 0")

    try:
        # find this atomic number's symbol, by looking at the list of atomic symbols
        return atomic_symbols[number - 1]
    except IndexError:
        # if this atomic number was out of range, then it is an invalid atomic number
        raise InvalidValueError("atomic number", number, "less than {}".format(len(atomic_symbols))) from None

def symbol_to_mass(symbol):
    """
    Finds the atomic mass of the atom with a given atomic symbol

    Args:
        symbol              - The 1 or 2 letter atomic symbol to find the mass of. For example: "He", "F".
                Case non-sensitive.

    Returns:
        The atomic mass of the given atom.
    """

    # the atomic mass is indexed in the list atomic_masses as the symbol's atomic number minus 1
    return atomic_masses[symbol_to_number(symbol) - 1]

def symbol_to_radius(symbol):
    """
    Finds the atomic radius of the atom with a given atomic symbol in angstroms.

    Args:
        symbol              - The 1 or 2 letter atomic symbol to find the radius of. For example: "He", "F".
                Case non-sensitive.

    Returns:
        The atomic radius of the given atom in angstroms.
    """

    # the atomic radius is indexed in the list atomic_radii as the symbol's atomic number minus 1
    return atomic_radii[symbol_to_number(symbol) - 1]

def symbol_to_covalent_radius(symbol):
    """
    Finds the covalent radius of the atom with a given atomic symbol in angstroms.

    Covalent radius is half the distance between 2 singly bonded atoms of the same type.

    Args:
        symbol              - The 1 or 2 letter atomic symbol to find the radius of. For example: "He", "F".
                Case non-sensitive.

    Returns:
        The covalent radius of the given atom in angstroms.
    """

    # the atomic covalent radius is indexed in the list covalent_radii as the symbol's atomic number minus 1
    return covalent_radii[symbol_to_number(symbol) - 1]


def symbol_to_vdw_radius(symbol):
    """
    Finds the van der Waals radius of the atom with a given atomic symbol in angstroms.

    Van der Waals radius is defined as half of the internuclear separation of two 
    non-bonded atoms of the same element on their closest possible approach.

    Args:
        symbol - The 1 or 2 letter atomic symbol to find the radius of. For example: "He", "F". Case non-sensitive.

    Returns:
        The van der Waals radius of the given atom in angstroms, if defined. Otherwise, throws an Exception. 
    """

    #the van der waals radius is indexed in the list vdw_radii as the symbol's atomic number minus 1
    vdw_radius = vdw_radii[symbol_to_number(symbol) - 1]

    if vdw_radius != -1.0:
        return vdw_radius
    else:
        raise InvalidValueError("Element", symbol, "has no valid van der Waals radius defined!") 


def symbol_to_free_polarizability(symbol):
    """
    Finds the free polarizability, taken from the published reference, of the atom with a given atomic symbol in angstrom^3.

    Args:
        symbol - The 1 or 2 letter atomic symbol to find the radius of. For example: "He", "F". Case non-sensitive.

    Returns:
        The free polarizability of the given atom in angstrom^3.
    """

    #the free polarizability is indexed in the list free_polarizabilities as the symbol's atomic number minus 1
    return bohr_to_ang**3 * free_polarizabilities[symbol_to_number(symbol) - 1]

def symbol_to_ccsdt_free_polarizability(symbol):
    """
    Finds the free polarizability, calculated using MolPRO with ccsd(t), of the atom with a given atomic symbol in angstrom^3.

    Args:
        symbol - The 1 or 2 letter atomic symbol to find the radius of. For example: "He", "F". Case non-sensitive.

    Returns:
        The free polarizability of the given atom in angstrom^3, if calculated. Otherwise, throws an Exception. 
    """

    try:
        return ccsdt_free_polarizabilities[symbol]
    except KeyError:
        raise InvalidValueError("Atom name", symbol, "free polarizability not calcualted for this atom; must be one of {}".format(", ".join(ccsdt_free_polarizabilities.keys())))

