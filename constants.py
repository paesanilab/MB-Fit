from exceptions import InvalidValueError

avogadro = 6.02214129E+23    # NIST 2010
au_to_ev = 27.21138505       # NIST 2010
au_to_joule = 4.35974434E-18 # NIST 2010
cal_to_joule = 4.184
bohr = 0.52917721067e-10 # meter, 2014 CODATA
bohr_to_ang = bohr * 1.0e10 # 0.52 ang / 1 bohr
ang_to_bohr = 1/bohr_to_ang # 1.82 bohr / 1 ang
# derived                                                                      
au_to_kcal = au_to_joule * avogadro * 1E-03 / cal_to_joule
ev_to_kcal = au_to_kcal / au_to_ev
au_per_bohr6_to_kcal_per_ang6 = au_to_kcal * (bohr_to_ang ** 6) # I guess this is right?

# atomic symbols ordered from lowest atomic number to highest
atomic_symbols = [
    "H",                                                                                                                                    "He",
    "Li",   "Be",                                                                                   "B",    "C",    "N",    "O",    "F",    "Ne",
    "Na",   "Mg",                                                                                   "Al",   "Si",   "P",    "S",    "Cl",   "Ar",
    "K",    "Ca",   "Sc",   "Ti",   "V",    "Cr",   "Mn",   "Fe",   "Co",   "Ni",   "Cu",   "Zn",   "Ga",   "Ge",   "As",   "Se",   "Br",   "Kr",
]

atomic_masses = [
    1.008,                                                                                                                                  4.0026,
    6.94,   9.0122,                                                                                 10.81,  12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305,                                                                                 26.982, 28.085, 30.974, 32.06,  35.45,  39.948,
    39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38,  69.723, 72.630, 74.922, 78.971, 79.904, 83.798,
]

def symbol_to_number(symbol):
    """
    Converts an atomic symbol to an atomic number
    """

    # Make first letter uppercase and second letter (if exists) lowercase
    symbol = symbol[:1].upper() + symbol[1:].lower()

    try:
        return atomic_symbols.index(symbol) + 1
    except ValueError:
        raise InvalidValueError("atomic symbol", symbol, "a valid 1 or 2 letter atomic symbol") from None

def number_to_symbol(number):
    try:
        if number < 1:
            raise ValueError("{} is not a recognized atomic number".format(number))
        return atomic_symbols[number - 1]
    except IndexError:
        raise InvalidValueError("atomic number", number, "less than {}".format(len(atomic_symbols))) from None

def symbol_to_mass(symbol):

    return atomic_masses[symbol_to_number(symbol) - 1]
