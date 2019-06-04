from potential_fitting.exceptions import XYZFormatError, InconsistentValueError
from .molecule import Molecule

'''
Generates a list of Molecule objects from xyz files in the given directory
'''
def xyz_to_molecules(file_path, settings = None):

    if settings is None:
        charge_per_fragment = [0]
        spin_per_fragment = [1]
        name_per_fragment = ["noname"]
        
        with open(file_path, "r") as xyz_file:
            total_atoms = int(xyz_file.readline())

        atoms_per_fragment = [total_atoms]        

        symmetry = ""

        symmetry_class = 65

        # loop over each atom assigning it a unique symmetry class
        for atom_index in range(total_atoms):

            symmetry += "{}1".format(chr(symmetry_class))

            symmetry_class += 1


        symmetry_per_fragment = [symmetry]

        SMILE = ""

        with open(file_path, "r") as xyz_file:
            for line in xyz_file.readlines()[2:2 + total_atoms]:
                SMILE += line.split()[0]

        SMILE_per_fragment = [SMILE]

    else:
        # the for loop just changes the string array to an int array
        atoms_per_fragment = [int(atom_count) for atom_count in settings.get("molecule", "fragments").split(",")]
        charge_per_fragment = [int(charge) for charge in settings.get("molecule", "charges").split(",")]
        spin_per_fragment = [int(spin) for spin in settings.get("molecule", "spins").split(",")]
        name_per_fragment = settings.get("molecule", "names").split(",")
        symmetry_per_fragment = settings.get("molecule", "symmetry").split(",")
        SMILE_per_fragment = settings.get("molecule", "SMILES").split(",")
        

    # open the file
    with open(file_path, "r") as xyz_file:

        # define the list of molecules
        molecules = []

        while True:
            try:
                molecule = Molecule.read_xyz_file(xyz_file, atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_per_fragment, symmetry_per_fragment, SMILE_per_fragment)
                molecules.append(molecule)
            except StopIteration:
                break
        
        return molecules

def parse_training_set_file(file_path, settings = None):

    if settings is None:
        charge_per_fragment = [0]
        spin_per_fragment = [1]
        name_per_fragment = ["noname"]
        
        with open(file_path, "r") as xyz_file:
            total_atoms = int(xyz_file.readline())

        atoms_per_fragment = [total_atoms]        

        symmetry = ""

        symmetry_class = 65

        # loop over each atom assigning it a unique symmetry class
        for atom_index in range(total_atoms):

            symmetry += "{}1".format(chr(symmetry_class))

            symmetry_class += 1


        symmetry_per_fragment = [symmetry]

    else:
        # the for loop just changes the string array to an int array
        atoms_per_fragment = [int(atom_count) for atom_count in settings.get("molecule", "fragments").split(",")]
        charge_per_fragment = [int(charge) for charge in settings.get("molecule", "charges").split(",")]
        spin_per_fragment = [int(spin) for spin in settings.get("molecule", "spins").split(",")]
        name_per_fragment = settings.get("molecule", "names").split(",")
        symmetry_per_fragment = settings.get("molecule", "symmetry").split(",")
        

    # open the file
    with open(file_path, "r") as xyz_file:

        while True:
            try:
                yield Molecule().read_xyz_file(xyz_file, atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_per_fragment, symmetry_per_fragment)
                
            except StopIteration:
                break