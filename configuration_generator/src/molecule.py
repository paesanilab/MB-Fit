# molecule.py
#
# definition for the Molecule class that holds data for molecules in nmcgen

import numpy as np

class Molecule:    
    def __init__(self, xyz, au_conversion=1.0):
        self.read_xyz(xyz, au_conversion)
    
    def read_xyz(self, xyz, au_conversion=1.0):
        self.num_atoms = 0
        atom_str = ""
        
        for i, line in enumerate(xyz.splitlines()):
            if i > 1:
                atom_str += line + "\n"
                if line.strip():
                    self.num_atoms += 1
        
        self.atoms = [""] * self.num_atoms
        self.coordinates = np.zeros((self.num_atoms, 3))
        
        for i, line in enumerate(atom_str.splitlines()):
            if line.strip():
                parts = line.split()
                
                atom = parts[0]
                
                if len(atom) > 1:
                    atom = atom[:1] + atom[-1:].lower()
                    
                self.atoms[i] = atom
                self.coordinates[i] = [float(parts[1])/au_conversion, float(parts[2])/au_conversion, float(parts[3])/au_conversion]
                
    def __str__(self):
        result = ""
        
        for i in range(self.num_atoms):
            result += self.atoms[i] + " "
            result += str(self.coordinates[i][0]) + " "
            result += str(self.coordinates[i][1]) + " "
            result += str(self.coordinates[i][2]) + "\n"
                
        return result
