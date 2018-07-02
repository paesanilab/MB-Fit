# molecule.py
#
# definition for the Molecule class that holds data for molecules in nmcgen

import psi4

class Molecule:    
    def __init__(self, charge, multiplicity, xyz):
        self.charge = charge
        self.multiplicity = multiplicity
        
        self.read_xyz(xyz)
    
    def read_xyz(self, xyz):
        self.atoms = ""
        self.num_atoms = 0
        
        for i, line in enumerate(xyz.splitlines()):
            if i > 1:
                self.atoms += line + "\n"
                if line.strip():
                    self.num_atoms += 1
                
    def read_psi4_mol(self, psi4_mol):
        self.read_xyz(psi4_mol.create_psi4_string_from_molecule())
    
    def psi4(self):
        geo = self.atoms + "\n" + str(self.charge) + " " + str(self.multiplicity)
        return psi4.geometry(geo)
