'''
The molecule class, which holds info for the molecule and the atoms
'''
class molecule:

    '''
    Ctor: Create this molecule from an array of atoms
    '''
    def __init__(self, atomArr):
        self.atoms = []
        for atom in atomArr:
            self.atoms.append(atom)

    '''
    Getter for list of atoms
    '''
    def getAtoms(self):
        return self.atoms

    '''
    Getter for molecule size (number of atoms)
    '''
    def getSize(self):
        return len(self.atoms)

    '''
    Setter for single-point energy of molecule
    '''
    def setEnergy(self, energy):
        self.energy = energy

    '''
    Getter for single-point energy of molecule
    '''
    def getEnergy(self):
        return self.energy

    '''
    String representation of the molecule
    '''
    def toString(self):

        # Concatenate the atom toString outputs
        atomStr = ""
        for atom in self.atoms:
            atomStr += atom.toString() + "\n"
        return atomStr

'''
The atom class, which holds info for atoms
'''
class atom:

    '''
    Ctor #1: Create this atom from an array of info
    '''
    def __init__(self, symbolArr):
        self.symbol = symbolArr[0]
        self.x = float(symbolArr[1])
        self.y = float(symbolArr[2])
        self.z = float(symbolArr[3])
    
    '''
    Ctor #2: Create this atom from manual input
    def __init__(self, symbol, x, y, z):
        self.symbol = symbol
        self.x = x
        self.y = y
        self.z = z
    '''

    '''
    Getter for symbol
    '''
    def getSymbol(self):
        return symbol

    '''
    Setter for symbol
    '''
    def setSymbol(self, newSym):
        self.symbol = newSym

    '''
    Getter for x
    '''
    def getX(self):
        return self.x

    '''
    Getter for y
    '''
    def getY(self):
        return self.y

    '''
    Getter for z
    '''
    def getZ(self):
        return self.z

    '''
    Setter for x
    '''
    def setX(self, newX):
        self.x = newX

    '''
    Setter for y
    '''
    def setY(self, newY):
        self.y = newY

    '''
    Setter for z
    '''
    def setZ(self, newZ):
        self.z = newZ

    '''
    Return a string representation of the atom
    '''
    def toString(self):
        return "{} {} {} {}".format(self.symbol, self.x, self.y, self.z)
