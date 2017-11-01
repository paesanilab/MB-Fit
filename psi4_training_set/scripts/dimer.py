'''
Extension class that deals with dimers
'''
class dimer:

    def __init__(mol1, mol2):
        self.mol1 = mol1
        self.mol2 = mol2

    def setMol1(mol1):
        self.mol1 = mol1

    def setMol2(mol2):
        self.mol2 = mol2

    def getMol1():
        return self.mol1

    def getMol2():
        return self.mol2

    def setEnergy(energy):
        self.energy = energy

    def getEnergy(energy):
        return self.energy

    '''
    Interaction energy!
    '''
    def intEnergy():
        return self.energy - self.mol1.getEnergy() - self.mol2.getEnergy()
