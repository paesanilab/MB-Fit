class PeriodicTable(object):
    def __init__(self):
        self.atoms = []

    def get_number(self, number):
        return self.atoms[number]

    def get_symbol(self, symbol):
        for atom in self.atoms:
            if atom.symbol == symbol:
                return atom
    
        return None # should probably raise an exception instead

class Atom(object):
    def __init__(self, number, symbol):
        self.number = number
        self.symbol = symbol
        self.radius = None
        self.covalent_radius = None
        self.mass = {}

    def add_isotope(self, number, mass, weight):
        self.mass[number] = (mass, weight)

    def get_average_mass(self):
        weighted_mass = 0
        total_weight = 0
        for mass, weight in self.mass:
            weighted_mass += mass * weight
            total_weight += weight

        return weighted_mass / total_weight

periodic_table = PeriodicTable()
