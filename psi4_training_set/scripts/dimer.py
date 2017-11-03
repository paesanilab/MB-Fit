class Dimer:
    """ Extension class that deals with dimers """

    def __init__(self, molecules):
        """Dimer constructor"""
        self.mol1 = molecules[0]
        self.mol2 = molecules[1]

    def set_mol1(self, mol1):
        """ Sets the first molecule within the dimer """
        self.mol1 = mol1

    def set_mol2(self, mol2):
        """ Sets the second molecule within the dimer """
        self.mol2 = mol2

    def set_energy(self, energy):
        """ Sets the energy of this dimer """
        self.energy = energy

    def __repr__(self):
        """ String representation of dimer """
        return str(self.mol1) + str(self.mol2)

    def int_energy():
        """ Interaction energy fo a dimer.
            This is defined by E(AB) - E(A) - E(B), 
            A and B representing molecules
        """
        return self.energy - self.mol1.energy - self.mol2.energy
