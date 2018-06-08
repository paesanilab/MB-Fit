from molecule import Atom, Fragment, Molecule
from exceptions import InvalidFormatException
import configparser
import os

'''
Generates a list of Molecule objects from xyz files in the given directory
'''
def xyz_to_molecules(f, config):

    # Keep in mind that you're passing in a file from initializer, not a
    # directory!

    # number of atoms in each fragment
    atoms_per_fragment = []
    # charge of each fragment
    charge_per_fragment = []
    # spin of each fragment
    spin_per_fragment = []

    # get list of all files and subdirectories in directory
    #dir_contents = os.listdir(directory)

    # Don't try to find the settings.ini, we should now REQUIRE it instead
    '''
    if "settings.ini" in dir_contents:
        # create config object
        config = configparser.SafeConfigParser(allow_no_value=False)

        # read the conflig file into the configuration
        config.read(directory + "/settings.ini")
    '''
        
    # the for loop just changes the string array to an int array
    atoms_per_fragment = [int(atom_count) for atom_count in config["molecule"]["fragments"].split(",")]
    charge_per_fragment = [int(charge) for charge in config["molecule"]["charges"].split(",")]
    spin_per_fragment = [int(spin) for spin in config["molecule"]["spins"].split(",")]

    # define the list of molecules
    molecules = []

    # Duplicate file openings
    '''
    # loop thru all files in directory
    for file_name in dir_contents:
        # check if file_name is an xyz file
        if file_name[-4:] == ".xyz":
            # parse the xyz file
            xyz = open(directory + "/" + file_name, "r")
    '''
    while True:

        # read first line of file, should be num of atoms
        atom_num_line = readline(f)
        # check for end of input
        if atom_num_line == "":
            break

        # check for consistance between number of atoms in xyz file
        # and number of atoms in fragments array
        if int(atom_num_line) != sum(atoms_per_fragment):
            raise Exception("Atom count specified in xyz file does not match count in ini file")
        # read and discard the comment line
        readline(f, True)

        # init Molecule object
        molecule = Molecule()

        # start reading atom lines
        for atom_count, charge, spin in zip(atoms_per_fragment, charge_per_fragment, spin_per_fragment):

            # init Fragment object
            fragment = Fragment(charge, spin)

            for i in range(0, atom_count):
                print(i)
                # read a line
                atom_line = readline(f, True)
                
                # extract atom info from line
                symbol, x, y, z = atom_line.split()

                # add atom to fragment
                fragment.add_atom(Atom(symbol, float(x), float(y), float(z)))
            
            # add fragment to molecule
            molecule.add_fragment(fragment)

        # add molecule to list of molecules
        molecules.append(molecule)

    return molecules

def readline(xyz, error = False):
    value = xyz.readline()
    print(value)
    if value == "":
        if error:
            raise InvalidFormatException(xyz.name, "<PLACEHOLDER>", "Parsing ended in the middle of a molecule.")
    return value

'''
molecules = xyz_to_molecules("config_files")
for molecule in molecules:
    # why is there an extra newline?
    print("MOLECULE: \n")
    print(molecule.to_xyz([0, 1, 2]))
'''
