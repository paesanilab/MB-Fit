'''
Extension class that deals with dimers
'''
class dimer:

    '''
    Dimer constructor
    '''
    def __init__(self, molecules):
        self.mol1 = molecules[0]
        self.mol2 = molecules[1]

    '''
    Sets molecule #1
    '''
    def setMol1(self, mol1):
        self.mol1 = mol1

    '''
    Sets molecule #2
    '''
    def setMol2(self, mol2):
        self.mol2 = mol2

    '''
    Gets molecule #1
    '''
    def getMol1(self):
        return self.mol1

    '''
    Gets molecule #2
    '''
    def getMol2(self):
        return self.mol2

    '''
    Sets the energy of this dimer
    '''
    def setEnergy(self, energy):
        self.energy = energy

    '''
    Gets the energy of this dimer
    '''
    def getEnergy(self):
        return self.energy

    '''
    String representation of dimer
    '''
    def toString(self):
        return self.mol1.toString() + "\n" + self.mol2.toString()

    '''
    Interaction energy!
    '''
    def intEnergy():
        return self.energy - self.mol1.getEnergy() - self.mol2.getEnergy()
