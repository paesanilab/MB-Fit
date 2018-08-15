import os, sys
from molecule import Atom, Fragment, Molecule

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
from exceptions import XYZFormatError, InconsistentValueError

'''
Generates a list of Molecule objects from xyz files in the given directory
'''
def xyz_to_molecules(file_path, settings):

    # open the file
    with open(file_path, "r") as xyz_file:
        # the for loop just changes the string array to an int array
        atoms_per_fragment = [int(atom_count) for atom_count in settings.get("molecule", "fragments").split(",")]
        charge_per_fragment = [int(charge) for charge in settings.get("molecule", "charges").split(",")]
        spin_per_fragment = [int(spin) for spin in settings.get("molecule", "spins").split(",")]
        name_per_fragment = settings.get("molecule", "names").split(",")

        # define the list of molecules
        molecules = []

        while True:
            molecule = Molecule().read_xyz(xyz_file, name_per_fragment, atoms_per_fragment, charge_per_fragment, spin_per_fragment)
            if molecule is None: # need some way of knowing if the file is exhuasted, is raising an error inside .read_xyz() better?
                break
            molecules.append(molecule)
        
        return molecules
