import sys, os

from hashlib import sha1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
from exceptions import XYZFormatError, InvalidValueError, InconsistentValueError

class Atom(object):
    """
    A class for an atom in a fragment, stores it's charge, number of unpaired
    electrons, and coordinates
    """

    '''
    Initialize a new atom from the given information
    '''
    def __init__(self, name, x, y, z):
        # Single letter symbol corresponding to periodic table abbrehviation.
        # For example: H, O, N, Cl, He
        self.name = name
        # x position in angstroms
        self.x = x
        # y position in angstroms
        self.y = y
        # z position in angstoms
        self.z = z

    '''
    Get the name of this symbol as its periodic table abbreviation
    '''
    def get_name(self):
        return self.name

    '''
    Get the x of this atom
    '''
    def get_x(self):
        return self.x

    '''
    Get the y of this atom
    '''
    def get_y(self):
        return self.y

    '''
    Get the z of this atom
    '''
    def get_z(self):
        return self.z

    '''
    Returns a string representing the information in this atom in the xyz file
    format
    '''
    def to_xyz(self):
        return "{:2} {:22.14e} {:22.14e} {:22.14e}".format(self.name, self.x, self.y, self.z)

    '''
    Returns a string representing the information in this atom in the xyz file
    format, specifying this atom as a ghost atom
    '''
    def to_ghost_xyz(self):
        return "@{:2} {:22.14e} {:22.14e} {:22.14e}".format(self.name, self.x, self.y, self.z)

class Fragment(object):
    """
    A class for a fragment of a Molecule. Contains atoms
    """
    
    '''
    Initialize a new Fragment with an empty atoms list
    '''
    def __init__(self, name, charge, spin_multiplicity):
        self.name = name
        # Array of atoms in this molecule
        self.atoms = []
        # charge of this fragment
        self.charge = charge
        # spin_multiplicity multiplicity of this fragment
        if spin_multiplicity < 1:
            raise InvalidValueError("spin multiplicity", spin_multiplicity, "1 or greater")
        self.spin_multiplicity = spin_multiplicity

    '''
    Gets the name of this fragment
    '''
    def get_name(self):
        return self.name

    '''
    Add an Atom to this Fragment.
    '''
    def add_atom(self, atom):
        self.atoms.append(atom)

    '''
    Gets a list of the atoms in this Fragment, sorted in standard order
    '''
    def get_atoms(self):
        return sorted(self.atoms, key=lambda atom: atom.to_xyz())

    '''
    Get the total charge of this fragment.
    '''
    def get_charge(self):
        return self.charge

    '''
    Get the spin_multiplicity multiplicity of this fragment
    '''
    def get_spin_multiplicity(self):
        return self.spin_multiplicity

    '''
    Gets the number of atoms in this fragment
    '''
    def get_num_atoms(self):
        return len(self.atoms)
 
    '''
    Returns a string representing the information in this fragment in the xyz
    file format in standard order
    '''
    def to_xyz(self):
        """ 
        Builds a string used for an output.
        """
        string = ""
        for atom in self.get_atoms():
            string += atom.to_xyz() + "\n"
        return string

    '''
    Returns a string represeting the information in this fragment in the xyz
    file format specifying the atoms in this fragment as ghost atoms in standard order
    '''
    def to_ghost_xyz(self):
        string = ""
        for atom in self.get_atoms():
            string += atom.to_ghost_xyz() + "\n"
        return string
    
    '''
    Initializes this fragment from a string, the returns itself
    '''
    def read_xyz(self, string):
        # split the string into an array of lines
        lines = string.splitlines(True)

        for line in lines:
            # construct each atom from the info in this line
            try:
                symbol, x, y, z = line.split()
            except ValueError:
                raise XYZFormatError(line, "ATOMIC_SYMBOL X Y Z") from None
            self.add_atom(Atom(symbol, float(x), float(y), float(z)))

        return self

class Molecule(object):
    """
    A Molecule holds an array of fragments.
    """

    '''
    Construct a new Molecule.
    '''
    def __init__(self):
        # list of fragments in this molecule
        self.fragments = []
        # list of energies for this molecule, filled in by get_nmer_energies
        self.energies = {}
        # list of nmer_energies for this molecule, filled by get_nmer_energies
        self.nmer_energies = []

        self.mb_energies = []

    '''
    gets the name of this molecule
    '''
    def get_name(self):
        return "-".join([fragment.get_name() for fragment in self.get_fragments()])

    '''
    Add a Fragment to this Molecule.
    '''
    def add_fragment(self, fragment):
        self.fragments.append(fragment)

    '''
    Gets a list of all the fragments in this molecule in standard order
    '''
    def get_fragments(self):
        return sorted(self.fragments, key=lambda fragment: fragment.get_name() + " " + fragment.to_xyz())

    '''
    Gets a list of all the atoms in this molecule in standard order
    '''
    def get_atoms(self):
        atoms = []
        for fragment in self.get_fragments():
            atoms += fragment.get_atoms()
        return atoms

    '''
    Get total charge of this Molecule by summing charges of Fragments.
    '''
    def get_charge(self, fragments = None):
        if fragments == None:
            fragments = range(len(self.get_fragments()))
        charge = 0
        for index in fragments:
            charge += self.get_fragments()[index].get_charge()
        return charge

    '''
    Get total spin_multiplicity multiplicity of this Molecule by summing spin_multiplicity
    of Fragments.
    '''
    def get_spin_multiplicity(self, fragments = None):
        if fragments == None:
            fragments = range(len(self.get_fragments()))
        spin_multiplicity = 1
        for index in fragments:
            spin_multiplicity += self.get_fragments()[index].get_spin_multiplicity() - 1
        return spin_multiplicity

    '''
    Gets the number of Fragments in this Molecule
    '''
    def get_num_fragments(self):
        return len(self.get_fragments())

    '''
    Gets the number of Atoms in this Molecule
    '''
    def get_num_atoms(self):
        atoms = 0
        for fragment in self.get_fragments():
            atoms += fragment.get_num_atoms()
        return atoms
    

    '''
    Returns a string representing the fragments of this Molecule in standard order

    Fragments should be specified by STANDARD ORDER
    '''
    def to_xyz(self, fragments = None, cp = False):
        # by default, use all fragments
        if fragments == None:
            fragments = range(self.get_num_fragments())
        string = ""
        for index in range(len(self.get_fragments())):
            if index in fragments:
                string += self.get_fragments()[index].to_xyz()
            elif cp:
                string += self.get_fragments()[index].to_ghost_xyz() 
        return string[:-1] # removes last character of string (extra newline)

    '''
    Returns a string containing indicies and energies of nbody fragment
    combinations in the format of the log file
    '''
    def log_frag_energy(self):
        string = ""
        # for each item in energies, add its combination indicies and energy
        # to the output string
        for combination in self.energies.keys():
            string += "E{}: {}\n".format(combination, "%.8f"%self.energies[combination])
        return string

    '''
    Returns a string containing the many body interaction energies, in the
    format of the log file
    '''
    def log_mb_energy(self, limit):
        string = ""
        
        for index in range(limit):
            string += "V_{}B: {}\n".format(index + 1, "%.8f"%self.mb_energies[index])
        
        return string
    

    '''
    Clears the energies, nmer_energies, and mb_energies fields to make way for
    new calculations
    '''
    def clear(self):
        self.energies = {}
        self.nmer_energies = []
        self.mb_energies = []

    '''
    Generates a SHA1 hash based on our molecule.
    We are using symbols, coordinates, and charges.
    Spin multiplicity to be added later.
    '''
    def get_SHA1(self):
        
        hash_string = self.to_xyz() + "\n" + str(self.get_charge()) + "\n" + str(self.get_spin_multiplicity())
        return sha1(hash_string.encode()).hexdigest()

    '''
    Initializes this molecule from a string, the returns itself
    '''
    def read_xyz(self, string, names, atoms_per_fragment, charges, spin_multiplicities):
        if not len(atoms_per_fragment) == len(charges):
            raise InconsistentValueError("atoms per fragment", "charges per fragment", atoms_per_fragment, charges, "lists must be same length")
        if not len(atoms_per_fragment) == len(spin_multiplicities):
            raise InconsistentValueError("atoms per fragment", "spin multiplicities per fragment", atoms_per_fragment, spin_multiplicities, "lists must be same length")
        if not len(atoms_per_fragment) == len(names):
            raise InconsistentValueError("atoms per fragment", "fragment names", atoms_per_fragment, names, "lists must be same length")

        # split the string into an array of lines
        lines = string.splitlines(True)

        if len(lines) != sum(atoms_per_fragment):
            raise InconsistentValueError("atoms per fragment", "number of atoms in input xyz", atoms_per_fragment, len(lines), "number of atoms in input xyz should equal the sum of the terms of atoms per fragment")

        for atom_count, charge, spin_multiplicity in zip(atoms_per_fragment, charges, spin_multiplicities):
            # construct each fragment from its charge, spin multiplicity and its line's from the string
            self.add_fragment(Fragment(charge, spin_multiplicity).read_xyz("".join(lines[:atom_count])))
            lines = lines[atom_count:]

        return self
