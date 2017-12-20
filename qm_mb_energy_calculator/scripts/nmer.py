class Nmer:
    """ A class that attempts to make polymers as generic as possible. """

    def __init__(self, molecules):
        self.mols = molecules

    def size(self):
        """ Gives the total amount of molecules """
        return len(self.mols)

    def total_atoms(self):
        """ Gives the total atom amount of the dimer """
        total = 0
        for mol in self.mols:
            total += mol.total_atoms()
        return total

    def set_energy(self, energy):
        self.energy = energy

    def __repr__(self):
        ret_str = ""
        for mol in self.mols:
            ret_str += str(mol)
        ret_str = ret_str[:-1]
        return ret_str

    def nmer_comb_str(self, comb):
        """ Builds a string representation of a N-mer part based on an input
            array.
        """
        ret_str = ""
        for mol_index in comb:
            ret_str += str(self.mols[mol_index])
        return ret_str

    def energy_str(self):
        """ String representation of energies used in output """
        ret_str = "%.8f"%self.energy
        for mol in self.mols:
            ret_str += " " + mol.energy_str()
        return ret_str

