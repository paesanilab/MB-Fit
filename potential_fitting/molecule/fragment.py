import numpy
import itertools
import collections

from .atom import Atom

from potential_fitting.exceptions import InvalidValueError, InconsistentValueError, XYZFormatError

# absolute package imports
from potential_fitting.polynomials.molecule_in_parser import FragmentParser

class Fragment(object):
    """
    Stores name, charge, spin multiplicity, and atoms of a fragment
    """
    
    def __init__(self, name, charge, spin_multiplicity):
        """
        Creates a new Fragment

        Args:
            name - the name of this fragment ('H2O', 'CO2', etc)
            charge - the charge of this fragment
            spin_multiplicity - the spin mulitplicity of this fragment, should be greater than or equal to 1

        Returns:
            A new Fragment
        """

        self.name = name
        # Array of atoms in this molecule
        self.atoms = []
        # charge of this fragment
        self.charge = charge
        # spin_multiplicity multiplicity of this fragment
        if spin_multiplicity < 1:
            raise InvalidValueError("spin multiplicity", spin_multiplicity, "1 or greater")
        self.spin_multiplicity = spin_multiplicity

        self.index = -1

    def get_name(self):
        """
        Gets the name of this fragment

        Args:
            None

        Returns:
            The name of this fragment
        """

        return self.name

    def add_atom(self, atom):
        """
        Adds an atom to this fragment

        Args:
            atom    - the atom to add

        Returns:
            None
        """

        for existing_atom in self.get_atoms():
            if existing_atom.get_symmetry_class() == atom.get_symmetry_class() and existing_atom.get_name() != atom.get_name():
                raise InconsistentValueError("symmetry class of {}".format(existing_atom.get_name()), "symmetry class of {}".format(atom.get_name()), existing_atom.get_symmetry_class(), atom.get_symmetry_class(), "atoms with a different atomic symbol in the same fragment cannot have the same symmetry class.")

        self.atoms.append(atom)

    def get_atoms(self):
        """
        Gets a list of the atoms in this fragment, sorted in standard order (alphabetically by xyz representation

        Args:
            None

        Returns:
            A list of the atoms in this fragment, sorted in standard order
        """

        return sorted(self.atoms, key=lambda atom: atom.get_name() + atom.get_symmetry_class())

    def get_index(self):
        """
        Gets the index of this fragment in the molecule

        Args:
            None

        Returns:
            The index of this framgnet
        """

        return self.index

    def get_symmetry(self):
        """
        Gets the symmetry of this fragment

        Args:
            None

        Returns:
            the symmetry of this fragment in A1B2 form
        """

        # used to build the symmetry string
        symmetry = self.get_atoms()[0].get_symmetry_class()

        # used to count how many atoms there are of the current symmetry
        symmetric_atom_count = 1

        for atom in self.get_atoms()[1:]:

            # if this atom has a different symmetry than the one before it
            if atom.get_symmetry_class() != symmetry[-1]:

                # record the number of atoms with the previous symmetry
                symmetry += str(symmetric_atom_count) + atom.get_symmetry_class()

                # reset the atom counter
                symmetric_atom_count = 1

            # if this atom has the same symmetry than the one before it
            else:

                # record this atom in the atom count
                symmetric_atom_count += 1

        # record the number of atoms with the last symmetry
        symmetry += str(symmetric_atom_count)

        return symmetry

    def get_charge(self):
        """
        Gets the charge of this fragment

        Args:
            None

        Returns:
            The charge of this fragment
        """

        return self.charge

    def get_spin_multiplicity(self):
        """
        Gets the spin multiplicity of this fragment

        Args:
            None

        Returns:
            The spin multiplicity of this fragment
        """

        return self.spin_multiplicity

    def get_num_atoms(self):
        """
        Gets the number of atoms in this fragment

        Args:
            None

        Returns:
            The number of atoms in this fragment
        """

        return len(self.atoms)

    def translate(self, x, y, z):
        """
        Translates all the atoms in this fragment by the given coordinates

        Args:
            x   - amount to translate along x axis
            y   - amount to translate along y axis
            z   - amount to translate along z axis

        Returns:
            None
        """

        for atom in self.get_atoms():
            atom.translate(x, y, z)

    def rotate(self, quaternion, origin_x = 0, origin_y = 0, origin_z = 0):
        """
        Rotates this Fragment using the rotation defined by the given Quaternion

        Args:
            quaternion - the Quaternion to rotate by
            origin_x - x position of the point to rotate around, default is 0
            origin_y - y position of the point to rotate around, default is 0
            origin_z - z position of the point to rotate around, default is 0

        Returns:
            None
        """

        for atom in self.get_atoms():
            atom.rotate(quaternion, origin_x, origin_y, origin_z)

    def get_connectivity_matrix(self):

    	# construct a matrix of size n by n where n is the number of atoms in this fragment
        # a value of 1 in row a and column b means that atom a and b are bonded
        connectivity_matrix = [[0 for k in range(self.get_num_atoms())] for i in range(self.get_num_atoms())]

        # loop over each pair of atoms
        for index1, atom1 in enumerate(self.get_atoms()):
            for index2, atom2 in enumerate(self.get_atoms()[index1 + 1:]):
                index2 += index1 + 1

                # if these atoms are bonded, set their values in the connectivity matrix to 1.
                if atom1.is_bonded(atom2):
                    connectivity_matrix[index1][index2] = 1
                    connectivity_matrix[index2][index1] = 1

    def get_excluded_pairs(self, max_exclusion = 3):
        """
        Gets the excluded pairs lists for this fragment

        Args:
            max_exclusion - get the excluded pairs up to 1x where x is max_exclusion, defualt is 3

        Returns:
            a tuple consisting of (excluded_12, excluded_13, ..., excluded_1x) lists
        """

        excluded_pairs = []

        connectivity_matrix = self.get_connectivity_matrix()

        # current matrix represents connectivity_matrix^x where x is the same as as in the excluded_1x pairs we are currently generating
        current_matrix = connectivity_matrix

        excluded_pairs_12 = set()

        # loop over each pair of atoms
        for index1, atom1 in enumerate(self.get_atoms()):
            for index2, atom2 in enumerate(self.get_atoms()[index1 + 1:]):
                index2 += index1 + 1

                # if the value in the current matrix is at least 1, then these atoms are 1 bond apart, and are added to the excluded_pairs_12 list.
                if current_matrix[index1][index2] > 0:
                    excluded_pairs_12.add((index1, index2))

        # add the excluded_pairs_12 to the list of all excluded pairs
        excluded_pairs.append(excluded_pairs_12)

        for i in range(max_exclusion - 1):

            # current matrix is multiplied by connectivity_matrix so that each iteration of the loop current_matrix = connectivity_matrix^(i + 1)
            current_matrix = numpy.matmul(current_matrix, connectivity_matrix)

            excluded_pairs_1x = set()

            # loop over each pair of atoms
            for index1, atom1 in enumerate(self.get_atoms()):
                for index2, atom2 in enumerate(self.get_atoms()[index1 + 1:]):
                    index2 += index1 + 1

                    # if the value in the connectivity matrix is at least 1, then these atoms are x bonds apart, and are added to the excluded_pairs_1x list.
                    if current_matrix[index1][index2] > 0:
                        excluded_pairs_1x.add((index1, index2))

            # filter out all terms inside other excluded lists from the new excluded list
            for excluded_pairs_1y in excluded_pairs:
                excluded_pairs_1x -= excluded_pairs_1y

            # add the excluded_pairs_1x to the list of all excluded pairs
            excluded_pairs.append(excluded_pairs_1x)

        return [[list(pair) for pair in excluded_pairs_1x] for excluded_pairs_1x in excluded_pairs]


    def to_xyz(self):
        """ 
        Gets the string representation of this fragment in the xyz file format

        Args:
            None

        Returns:
            String containing the atomic symbols and coordinates of each atom in this fragment in standard order
        """

        string = ""

        # add each atom to the string
        for atom in self.get_atoms():
            string += atom.to_xyz() + "\n"

        return string

    def to_ghost_xyz(self):
        """ 
        Gets the string representation of this fragment in the xyz file format as a ghost fragment

        Args:
            None

        Returns:
            string containing the atomic symbols and coordinates of each atom in this fragment in standard order as ghost atoms
        """

        string = ""

        # add each atom to the string
        for atom in self.get_atoms():
            string += atom.to_ghost_xyz() + "\n"

        return string

    def read_xyz(self, string, symmetry):
        """
        Reads the given string into this Fragment with the given symmetry
        
        Args:
            string - the xyz file format string, just the lines with the atom information, no atom count / comment line
            symmetry - the symmetry of this Fragment

        Returns:
            self
        """

        # used to keep track of the index of the current atom symmetry class in the symmetry string
        class_index = 0
        # used to keep track of how many atoms of the current atom symmetry class have been made
        class_counter = 0

        fragment_parser = FragmentParser(symmetry, "a")

        lines = string.splitlines()

        for atom_type in fragment_parser.get_atom_types():

            for atom in atom_type.get_atoms():

                symmetry_class = "".join(itertools.takewhile(str.isupper, atom))

                try:
                    # read one line of the input string into an atom
                    symbol, x, y, z = lines[0].split()
                    lines = lines[1:]
                except ValueError:
                    # if we cannot parse the line, raise an error                
                    raise XYZFormatError(line, "ATOMIC_SYMBOL X Y Z") from None
                except IndexError:
                    # there are no more lines to read but are still symmetries
                    raise InconsistentValueError("atom lines in fragment xyz", "symmetry of fragment", len(string.splitlines()), symmetry, "sum of numbers in symmetry must equal number of atom lines in the fragment xyz")

                self.add_atom(Atom(symbol, symmetry_class, float(x), float(y), float(z)))

        # check if there are more lines than atoms in the symmetry
        if len(lines) != 0:
            raise InconsistentValueError("atom lines in fragment xyz", "symmetry of fragment", len(string.splitlines()), symmetry, "sum of numbers in symmetry must equal number of atom lines in the fragment xyz")

        return self

    def get_priority(self, connectivity_matrix, atom):

    priority = atom.get_mass()

    for

    return priority, next_atom


    def get_standard_order(self):

    	# make sure order of connectivity_matrix matches that of atoms list
    	connectivity_matrix = get_connectivity_matrix()

    	atoms = self.get_atoms()

    	atoms = sorted(atoms, lambda x: x.get_mass())

    	atom_names = [atom.get_name() for atom in atoms]

    	atom_name_dict = dict(collections.Counter(atom_names))

    	symmetry = ""
    	standard_order = [];

    	next_sym = "A"

    	for atom_name, count in atom_names_dict.items():
    		if count == 1:
    			symmetry += next_sym
    			next_sym = chr(ord(next_sym) + 1)
    			standard_order.append(atoms[0])
    			atoms = atoms[1:]
    		else:
    			tie_atoms = atoms[:count]
    			atoms = atoms[counts:]




    def __eq__(self, other):
        if not self.get_name() == other.get_name():
            return False
        if not self.get_charge() == other.get_charge():
            return False
        if not self.get_spin_multiplicity() == other.get_spin_multiplicity():
            return False

        for self_atom, other_atom in zip(self.get_atoms(), other.get_atoms()):
            if self_atom != other_atom:
                return False

        return True

    def __ne__(self, other):
        return not self == other
