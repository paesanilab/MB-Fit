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
        symmetry_per_fragment = settings.get("molecule", "symmetry").split(",")

        # define the list of molecules
        molecules = []

        while True:
            try:
                molecule = Molecule().read_xyz_file(xyz_file, atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_per_fragment, symmetry_per_fragment)
                molecules.append(molecule)
            except StopIteration:
                break
        
        return molecules
