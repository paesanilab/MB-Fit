class Molecule:

    def __init__(self, atom_arr):
        self.atoms = []
        for atom in atom_arr:
            self.atoms.append(atom)
    
    def size(self):
        return 1

    def total_atoms(self):
        return len(self.atoms)

    def set_energy(self, energy):
        """ Setter for single-point energy of molecule """
        self.energy = energy

    def nmer_comb_str(self, comb):
        return str(self)

    def __repr__(self):

        # Concatenate the atom toString outputs
        atom_str = ""
        for atom in self.atoms:
            atom_str += str(atom) + "\n"
        return atom_str

    def energy_str(self):
        return "%.8f"%self.energy

class Atom:

    def __init__(self, symbol_arr):
        """ Constructor that takes in an array of symbol and coordinates """
        self.symbol = symbol_arr[0]
        self.x = float(symbol_arr[1])
        self.y = float(symbol_arr[2])
        self.z = float(symbol_arr[3])
    
    def set_symbol(self, new_sym):
        self.symbol = new_sym

    def set_x(self, new_x):
        self.x = new_x

    def set_y(self, new_y):
        self.y = new_y

    def set_z(self, new_z):
        self.z = new_z

    def __repr__(self):
        return "{} {} {} {}".format(self.symbol, self.x, self.y, self.z)
