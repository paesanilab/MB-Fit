import os, sys
import configparser
from molecule import Atom, Fragment, Molecule

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
from exceptions import XYZFormatError, InconsistentValueError

'''
Generates a list of Molecule objects from xyz files in the given directory
'''
def xyz_to_molecules(f, config):

    # Keep in mind that you're passing in a file from initializer, not a
    # directory!

    # the for loop just changes the string array to an int array
    atoms_per_fragment = [int(atom_count) for atom_count in config["molecule"]["fragments"].split(",")]
    charge_per_fragment = [int(charge) for charge in config["molecule"]["charges"].split(",")]
    spin_per_fragment = [int(spin) for spin in config["molecule"]["spins"].split(",")]
    name_per_fragment = config["molecule"]["names"].split(",")

    if not len(atoms_per_fragment) == len(charge_per_fragment):
        raise InconsistentValueError("fragments", "charges", config["molecule"]["fragments"], config["molecule"]["charges"], "lists must be same length")
    if not len(atoms_per_fragment) == len(spin_per_fragment):
        raise InconsistentValueError("fragments", "spins", config["molecule"]["fragments"], config["molecule"]["spins"], "lists must be same length")
    if not len(atoms_per_fragment) == len(name_per_fragment):
        raise InconsistentValueError("fragments", "names", config["molecule"]["fragments"], config["molecule"]["names"], "lists must be same length")

    # define the list of molecules
    molecules = []

    while True:

        # read first line of file, should be num of atoms
        atom_num_line = f.readline()

        # check for end of input
        if atom_num_line == "":
            break

        # check for consistance between number of atoms in xyz file
        # and number of atoms in fragments array
        if int(atom_num_line) != sum(atoms_per_fragment):
            raise InconsistentValueError("total atoms in xyz file", "fragments", atom_num_line, config["molecule"]["fragments"], "fragments list must sum to total atoms from xyz file")

        # read and discard the comment line
        f.readline()

        # init Molecule object
        molecule = Molecule()

        # start reading atom lines
        for atom_count, name, charge, spin in zip(atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_per_fragment):

            # init Fragment object
            fragment = Fragment(name, charge, spin)

            for i in range(0, atom_count):
                # read a line
                atom_line = f.readline()
                if atom_line == "":
                    raise XYZFormatError("end of file reached while parsing", "expected additional atom line")
                
                # extract atom info from line
                try:
                    symbol, x, y, z = atom_line.split()
                except ValueError:
                    raise XYZFormatError(line, "ATOMIC_SYMBOL X Y Z") from None
                # add atom to fragment
                fragment.add_atom(Atom(symbol, float(x), float(y), float(z)))
            
            # add fragment to molecule
            molecule.add_fragment(fragment)

        # add molecule to list of molecules
        molecules.append(molecule)

    return molecules
