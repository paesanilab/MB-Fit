from molecule import Atom, Fragment, Molecule
from exceptions import InvalidFormatException

'''
Generates a list of Molecule objects from the given xyz input file
'''
def xyz_to_molecules(xyz):
    # define the list of Molecules
    molecules = []

    # first line of next Molecule to parse, defined outside because it needs to
    # be updated inside the loop
    firstline = readline(xyz)

    while True:
    
        # if firstline is empty, there are no more Molecules to parse
        if firstline == "":
            break
        # read total number of atoms from xyz file
        total_atoms = int(firstline)
        # read atoms per fragment from xyz file
        atoms_per_fragment = readline(xyz).split()

        # read atoms from xyz file
        atoms = []
        line = xyz.readline()
        while line and len(line.split()) == 4:
            atoms.append(line)
            line = xyz.readline()
        firstline = line
        

        # check to make sure total number of atoms, atoms per fragment, and given
        # atoms all indicate the same number of total atoms

        # sum the atoms that are supposed to be in each fragment
        fragment_atoms = 0
        for fragment in atoms_per_fragment:
            fragment_atoms +=   int(fragment)

        # check if the atoms in each fragment add up to the total number of atoms
        # from the first line of the xyz file
        if total_atoms != fragment_atoms:
            raise InvalidFormatException(xyz.name, "<PLACEHOLDER>", "Total atom number in input file ({}) does not match sum of fragment atoms ({}).".format(total_atoms, fragment_atoms))

        # check if the given number of atoms equals the total number of atoms from
        # the first line of the xyz file
        if total_atoms != len(atoms):
            raise InvalidFormatException(xyz.name, "<PLACEHOLDER>", "Total atom number in input file ({}) does not match actual number of atoms listed ({}).".format(total_atoms, len(atoms)))
        
        # if no errors were encountered, then input is valid and a Molecule object
        # can be generated

        # construct a new Molecule object
        molecule = Molecule()
        
        # loop through all the fragments, building them and adding them to the
        # Molecule
        for fragment_size in atoms_per_fragment:
            # construct a new fragment
            fragment = Fragment()
            # loop through each atom in this fragment, adding it to the fragment
            for index in range(int(fragment_size)):
                # String representation of the atom at the beginning of the atoms
                # list
                atom_str = atoms[0]
                atoms = atoms[1:]
                # split the atom_str into [name, x, y, z]
                atom_properties = atom_str.split()
                # name is first element of atom_properties
                name = atom_properties[0]
                # charge is 0 for now
                charge = 0
                # unpaired electrons is 0 for now
                unpaired_electrons = 0
                # coordinates are 2nd, 3rd, and 4th items of atom_properties
                x = float(atom_properties[1])
                y = float(atom_properties[2])
                z = float(atom_properties[3])
                # construct Atom from its properties
                atom = Atom(name, charge, unpaired_electrons, x, y, z)
                # add the Atom the the Fragment
                fragment.add_atom(atom)
            # add the Fragment to the Molecule
            molecule.add_fragment(fragment)
        # append the Molecule to the list of Molecules
        molecules.append(molecule)
    # return the completed list of molecules
    return molecules

def readline(xyz, error = False):
    value = xyz.readline()
    if value == "":
        if error:
            raise InvalidFormatException(xyz.name, "<PLACEHOLDER>", "Parsing ended in the middle of a molecule.")
    return value

'''
COMMENTED TEST CODE
molecules = xyz_to_molecules(open("trimer_input.xyz", 'r'))
for molecule in molecules:
    # why is there an extra newline?
    print("MOLECULE: \n")
    print(molecule.to_xyz([0, 1, 2]))
'''
