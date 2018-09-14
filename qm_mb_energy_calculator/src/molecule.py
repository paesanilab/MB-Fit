import sys, os
import numpy, math

from hashlib import sha1

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
from exceptions import XYZFormatError, InvalidValueError, InconsistentValueError
import constants

class Atom(object):
    """
    Stores name, x, y, and z of a single atom
    """

    def __init__(self, name, symmetry_class, x, y, z):
        """
        Creates a new atom

        Args:
            name    - the atomic symbol of this atom ('H', 'He', etc)
            x       - the x position of this atom in angstroms
            y       - the y position of this atom in angstroms
            z       - the z position of this atom in angstroms

        Returns:
            A new Atom
        """

        self.name = name
        self.symmetry_class = symmetry_class
        self.x = x
        self.y = y
        self.z = z

    def get_name(self):
        """
        Gets the name of this atom
        
        Args:
            None

        Returns:
            The atomic symbol of this atom ('H', 'He', etc)
        """

        return self.name

    def get_symmetry_class(self):
        """
        Gets the symmetry class of this atom

        Args:
            None

        Returns:
            The symmetry class of this atom
        """

        return self.symmetry_class

    def set_symmetry_class(self, symmetry_class):
        """
        Changes the symmetry class of this atom

        Args:
            symmetry_class - the new symmetry class of this atom

        Returns:
            None
        """

        self.symmetry_class = symmetry_class

    def get_number(self):
        """
        Gets the atomic number of this atom

        Args:
            None

        Returns:
            The atomic number of this atom
        """

        return constants.symbol_to_number(self.name)

    def get_mass(self):
        """
        Gets the atomic mass of this atom

        Args:
            None

        Returns:
            The atomic mass of this atom in g/mol
        """

        return constants.symbol_to_mass(self.name)

    def get_radius(self):
        """
        Gets the atomic radius of this atom

        Args:
            None

        Returns:
            The atomic radius of this atom
        """

        return constants.symbol_to_radius(self.name)

    def get_covalent_radius(self):
        """
        Gets the covalent radius of this atom

        Args:
            None

        Returns:
            The covalent radius of this atom
        """

        return constants.symbol_to_covalent_radius(self.name)
    def get_x(self):
        """
        Gets the x position of this atom

        Args:
            None

        Returns:
            The x position of this atom in angstroms
        """

        return self.x

    def get_y(self):
        """
        Gets the y position of this atom

        Args:
            None

        Returns:
            The y position of this atom in angstroms
        """

        return self.y

    def get_z(self):
        """
        Gets the z position of this atom

        Args:
            None

        Returns:
            The z position of this atom in angstroms
        """

        return self.z

    def set_x(self, x):
        """
        Sets the x position of this atom

        Args:
            x   - the new x position

        Returns:
            None
        """

        self.x = x

    def set_y(self, y):
        """
        Sets the y position of this atom

        Args:
            y   - the new y position

        Returns:
            None
        """

        self.y = y

    def set_z(self, z):
        """
        Sets the z position of this atom

        Args:
            z   - the new z position

        Returns:
            None
        """

        self.z = z

    def set_xyz(self, x, y, z):
        """
        Sets the x, y, and z positions of this atom

        Args:
            x   - the new x position
            y   - the new y position
            z   - the new z position

        Returns:
            None
        """

        self.x = x
        self.y = y
        self.z = z

    def translate(self, x, y, z):
        """
        Translates this atom by the given coordinates

        Args:
            x   - amount to translate along x axis
            y   - amount to translate along y axis
            z   - amount to translate along z axis

        Returns:
            None
        """

        self.x += x
        self.y += y
        self.z += z

    def rotate(self, quaternion, origin_x = 0, origin_y = 0, origin_z = 0):
        """
        Rotates this Atom using the rotation defined by the given Quaternion

        Args:
            quaternion - the Quaternion to rotate by
            origin_x - x position of the point to rotate around, default is 0
            origin_y - y position of the point to rotate around, default is 0
            origin_z - z position of the point to rotate around, default is 0

        Returns:
            None
        """

        x, y, z = quaternion.rotate(self.get_x(), self.get_y(), self.get_z(), origin_x, origin_y, origin_z)

        self.set_xyz(x, y, z)
        
    def distance(self, atom):
        """
        Finds the distance between another atom and this one

        Args:
            atom    - the atom to compare to this one

        Returns:
            distance between them, in Angstroms
        """
        # compute distance in 3d coordinate plane
        return math.sqrt((self.get_x() - atom.get_x()) ** 2 + (self.get_y() - atom.get_y()) ** 2 + (self.get_z() - atom.get_z()) ** 2)

    def to_xyz(self):
        """
        Gets the string representation of this atom in the xyz file format

        Args:
            None

        Returns:
            String containing this atom's atomic symbol and coordinates in the xyz format
        """
        return "{:2} {:22.14e} {:22.14e} {:22.14e}".format(self.name, self.x, self.y, self.z)

    def is_bonded(self, atom, bond_sensitivity = 1.1):
        """
        Calculates whether this atom is likely to be bonded to another based on their atomic radii and the distance between them.

        Args:
            atom    - the atom to check if this one is bonded to
            bond_sensitivity - the bond threshold is considered to be this * the sum of the atomic radii

        Returns:
            True if the distance between the atoms is less than bond_sensitivity * the sum of their covalent radii, otherwise False.
        """

        return self.distance(atom) < bond_sensitivity * (self.get_covalent_radius() + atom.get_covalent_radius())

    # NOTE: using @ prefix is not universal, setup some way to change ghost representation depending on platform.
    def to_ghost_xyz(self):
        """
        Gets the string representation of this atom in the xyz file format as a ghost atom

        Args:
            None

        Returns:
            String containing this atom's atomic symbol and coordinates in the xyz ghost atom format
        """
        return "@{}".format(self.to_xyz())

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

        return sorted(self.atoms, key=lambda atom: atom.get_symmetry_class())

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

    def get_excluded_pairs(self, max_exclusion = 3):
        """
        Gets the excluded pairs lists for this fragment

        Args:
            max_exclusion - get the excluded pairs up to 1x where x is max_exclusion, defualt is 3

        Returns:
            a tuple consisting of (excluded_12, excluded_13, ..., excluded_1x) lists
        """

        excluded_pairs = []

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

        return tuple(list(excluded_pairs_1x) for excluded_pairs_1x in excluded_pairs)


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
    
    def read_xyz(self, file, num_atoms, symmetry):
        """
        Reads a number of lines from an open file into this fragment

        Args:
            file  - the file to read lines from
            num_atoms - the number of lines (and atoms) to read from the xyz file

        Returns:
            This Fragment
        """

        symmetry_index = 0
        symmetry_count = 0


        # loop once for each atom expected to be read
        for i in range(num_atoms):

            if not symmetry_count < int(symmetry[2 * symmetry_index + 1]):
                symmetry_index += 1
                symmetry_count = 0

            # get the line corresponding to this atom
            line = file.readline()

            # if we have reached EOF, raise an error because we expected at least 1 more atom
            if line == "":
                raise XYZFormatError("end of file reached while parsing", "expected additional atom line")

            try:
                # construct each atom from the info in the line
                symbol, x, y, z = line.split()

            except ValueError:
                # if we cannot parse the line, raise an error                
                raise XYZFormatError(line, "ATOMIC_SYMBOL X Y Z") from None

            try:
                symmetry_class = symmetry[2 * symmetry_index]
            except IndexError:
                raise InconsistentValueError("num atoms in fragment", "symmetry of fragment", num_atoms, symmetry, "sum of numbers in symmetry must equal num atoms in fragment")

            self.add_atom(Atom(symbol, symmetry_class, float(x), float(y), float(z)))

            symmetry_count += 1

        return self

class Molecule(object):
    """
    Stores the fragments of a Molecule
    """

    def __init__(self):
        """
        Creates a new Molecule

        Args:
            None

        Returns:
            A new Molecule
        """

        # list of fragments in this molecule
        self.fragments = []
        # list of energies for this molecule, filled in by get_nmer_energies
        self.energies = {}
        # list of nmer_energies for this molecule, filled by get_nmer_energies
        self.nmer_energies = []

        self.mb_energies = []

    def get_name(self):
        """
        Gets the name of this molecule, consists of the names of the fragments in standard order connected by a dash '-'

        Args:
            None

        Returns:
            The name of this molecule
        """

        return "-".join([fragment.get_name() for fragment in self.get_fragments()])

    def get_symmetry(self):
        """
        Gets the symmetry of this molecule

        Args:
            None

        Returns:
            The symmetry of this molecule in A1B2_C1D1E1 form
        """

        # used to assemble the symmetry string
        symmetry = self.get_fragments()[0].get_symmetry()

        # add each fragment's symmetry to the string
        for fragment in self.get_fragments()[1:]:

            symmetry += "_" + fragment.get_symmetry()

        return symmetry

    def add_fragment(self, fragment):
        """
        Adds a fragment to this molecule

        Args:
            fragment    - the fragment to add

        Returns:
            None
        """

        # make sure the symmetry class of the atoms in this fragment doesn't violate the 1 symmetry class -> 1 atom type rule
        for existing_fragment in self.get_fragments():
            for atom1 in existing_fragment.get_atoms():
                for atom2 in fragment.get_atoms():
                    if atom1.get_symmetry_class() == atom2.get_symmetry_class() and atom1.get_name() != atom2.get_name():
                        raise InconsistentValueError("symmetry class of {}".format(atom1.get_name()), "symmetry class of {}".format(atom2.get_name()), atom1.get_symmetry_class(), atom2.get_symmetry_class(), "atoms with a different atomic symbol in the same molecule (even if different fragments) cannot have the same symmetry class.")

        self.fragments.append(fragment)

        # adjust the symmetry class of atoms to be in the right order
        next_symmetry = 65
        new_symmetries = {}
        for fragment in self.get_fragments():
            for atom in fragment.get_atoms():
                if atom.get_symmetry_class() not in new_symmetries:
                    new_symmetries[atom.get_symmetry_class()] = chr(next_symmetry)
                    next_symmetry += 1

        for fragment in self.get_fragments():
            for atom in fragment.get_atoms():
                atom.set_symmetry_class(new_symmetries[atom.get_symmetry_class()])

    def get_fragments(self):
        """
        Gets a list of the fragments in this molecule in standard order

        Args:
            None

        Returns:
            List of fragments in this molecule in standard order
        """

        return sorted(self.fragments, key=lambda fragment: fragment.get_name())

    def get_atoms(self):
        """
        Gets a list of the atoms in this molecule in standard order

        fragments are first sorted into standard order, and then atoms within those fragments are put in their standard order.

        Args:
            None

        Returns:
            List of atoms in this molecule in standard order
        """

        atoms = []
        for fragment in self.get_fragments():
            atoms += fragment.get_atoms()
        return atoms

    def get_charge(self, fragments = None):
        """
        Gets the charge of this molecule by summing the charges of its fragments

        Args:
            fragments   - list of fragment indicies; if specified, only get the charge of these fragments, default is to include all fragments

        Returns:
            Sum charge of all or some of the fragments of this molecule
        """

        if fragments == None:
            fragments = range(len(self.get_fragments()))
        charge = 0
        for index in fragments:
            charge += self.get_fragments()[index].get_charge()
        return charge

    def get_spin_multiplicity(self, fragments = None):
        """
        Gets the spin multiplicity of this molecule by summing the spin multiplicities of its fragments

        Args:
            fragments   - list of fragment indicies; if specified, only get the spin multiplicity of these fragments, default is to include all fragments

        Returns:
            Sum spin multiplicity of all or some of the fragments of this molecule
        """

        if fragments == None:
            fragments = range(len(self.get_fragments()))
        spin_multiplicity = 1
        for index in fragments:
            spin_multiplicity += self.get_fragments()[index].get_spin_multiplicity() - 1
        return spin_multiplicity

    def get_num_fragments(self):
        """
        Gets the number of fragments in this molecule

        Args:
            None

        Returns:
            Number of fragments in this molecule
        """

        return len(self.get_fragments())

    def get_num_atoms(self):
        """
        Gets the number of atoms in this molecule

        Args:
            None

        Returns:
            Number of atoms in this molecule
        """

        atoms = 0
        for fragment in self.get_fragments():
            atoms += fragment.get_num_atoms()
        return atoms

    def translate(self, x, y, z):
        """
        Translates all the atoms in this molecule by the given coordinates

        Args:
            x   - amount to translate along x axis
            y   - amount to translate along y axis
            z   - amount to translate along z axis
    
        Returns:
            None
        """

        for fragment in self.get_fragments():
            fragment.translate(x, y, z)

    def rotate(self, quaternion, origin_x = 0, origin_y = 0, origin_z = 0):
        """
        Rotates this Molecule using the rotation defined by the given Quaternion

        Args:
            quaternion - the Quaternion to rotate by
            origin_x - x position of the point to rotate around, default is 0
            origin_y - y position of the point to rotate around, default is 0
            origin_z - z position of the point to rotate around, default is 0

        Returns:
            None
        """

        for fragment in self.get_fragments():
            fragment.rotate(quaternion, origin_x, origin_y, origin_z)

    def move_to_center_of_mass(self):
        """
        Moves the molecule it its center of mass

        Args:
            None

        Returns:
            None
        """

        # keep track of the total weighted mass along each axis
        total_x = 0
        total_y = 0
        total_z = 0

        # keeps track of the total mass
        total_mass = 0

        # loop thru every atom in the molecule, adding its contribution to each coordinate mass
        for atom in self.get_atoms():
            total_x += atom.get_x() * atom.get_mass()
            total_y += atom.get_y() * atom.get_mass()
            total_z += atom.get_z() * atom.get_mass()

            total_mass += atom.get_mass()

        # calculate the center of mass my dividing the total weighted mass by the total mass
        center_x = total_x / total_mass
        center_y = total_y / total_mass
        center_z = total_z / total_mass

        # translate this molecule to the center of mass
        self.translate(-center_x, -center_y, -center_z)

    def rotate_on_principal_axes(self):
        """
        Rotates a molecule on to its principle axis

        Args:
            None

        Returns:
            None
        """

        # first we calculate the moment of inertia tensor
        # [ Ixx Ixy Ixz ]
        # [ Iyx Iyy Iyz ]
        # [ Izx Izy Izz ]
        I = [[0, 0, 0] for i in range(3)]

        # loop over every atom and add their contributions to the moment of inertia tensor
        for atom in self.get_atoms():
            # Ixx
            I[0][0] += (atom.get_y() ** 2 + atom.get_z() ** 2) * atom.get_mass()
            # Ixy
            I[1][0] += - (atom.get_x() * atom.get_y()) * atom.get_mass()
            # Ixz
            I[2][0] += - (atom.get_x() * atom.get_z()) * atom.get_mass()

            # Iyx
            I[0][1] += - (atom.get_y() * atom.get_x()) * atom.get_mass()
            # Iyy
            I[1][1] += (atom.get_x() ** 2 + atom.get_z() ** 2) * atom.get_mass()
            # Iyz
            I[2][1] += - (atom.get_y() * atom.get_z()) * atom.get_mass()

            # Izx
            I[0][2] += - (atom.get_z() * atom.get_x()) * atom.get_mass()
            # Izy
            I[1][2] += - (atom.get_z() * atom.get_y()) * atom.get_mass()
            # Izz
            I[2][2] += (atom.get_x() ** 2 + atom.get_y() ** 2) * atom.get_mass()

        # get numpy matrix from the matrix of principle moments
        inertia_tensor = numpy.matrix(I)

        # get the moments and principle axis as eigen values and eigen vectors
        (moments, principle_axes) = numpy.linalg.eig(inertia_tensor)

        # reorder the principle axes from largest eigen value to smallest
        idx = moments.argsort()[::-1]
        moments = moments[idx]
        principle_axes = principle_axes[:,idx]

        # update the position of each atom
        for atom in self.get_atoms():
            x, y, z = (numpy.matrix([atom.get_x(), atom.get_y(), atom.get_z()]) * principle_axes).getA1()
            atom.set_xyz(float(x), float(y), float(z))

    def rmsd(self, other):
        """
        Computes the RMSD between the positions of the atoms in two molecules

        molecules must have the same fragments and atoms or an InconsistentValueError will be raised.

        generally, you should make sure that both molecules have been moved to their center of mass and rotated on their principal axes.

        Args:
            other - the molecule to compare this one to

        Returns:
            The square-root of the mean squared distance between the atoms in this molecule and the other
        """

        # fist make sure these molecules have the same number of atoms
        if self.get_num_atoms() != other.get_num_atoms():
            raise InconsistentValueError("number of atoms in self", "number of atoms in other", self.get_num_atoms(), other.get_num_atoms(), "number of atoms in each molecule must be the same, make sure you are computing the rmsd of two molecules with the same atoms and fragments")

        squared_distance = 0

        # loop thru every pair of atoms in the two molecules
        for this_atom, other_atom in zip(self.get_atoms(), other.get_atoms()):

            # check to make sure that these atoms are the same type
            if this_atom.get_name() != other_atom.get_name():
                raise InconsistentValueError("self atom symbol", "other atom symbol", this_atom.get_name(), other_atom.get_name(), "symbols must be the same, make sure you are computing the rmsd of two molecules with the same atoms and fragments")

            # add this atom pair's contribution to the squared distance
            squared_distance += this_atom.distance(other_atom) ** 2

        # compute rmsd as sqrt of mean squared distance
        return math.sqrt(squared_distance / self.get_num_atoms())

    def distancermsd(self, other_molecule):
        """
        Computes the RMSD of intramolecular interatomic distances in the two molecules

        molecules must have the same fragments and atoms or an InconsistentValueError will be raised.

        generally, you should make sure that both molecules have been moved to their center of mass and rotated on their principal axes.

        Note:
            this function is distinct from rmsd() because this function takes the rmsd of the differneces between the distances between pairs of atoms within each molecule
            while rmsd() takes the rmsd of the distance between the positions of the same atoms in each molecule.

        Args:
            other_molecule - the molecule to ompare this one to

        Returns:
            the square-root of the mean squared difference in the distance between each pair of atoms in this molecule and the other
        """

        # fist make sure these molecules have the same number of atoms
        if self.get_num_atoms() != other_molecule.get_num_atoms():
            raise InconsistentValueError("number of atoms in self", "number of atoms in other", self.get_num_atoms(), other_molecule.get_num_atoms(), "number of atoms in each molecule must be the same, make sure you are computing the rmsd of two molecules with the same atoms and fragments")

        squared_distance_difference = 0

        # loop over each pair of atoms
        for atom_index, this_atom1, other_atom1 in zip(range(self.get_num_atoms()), self.get_atoms(), other_molecule.get_atoms()):
            for this_atom2, other_atom2 in zip(self.get_atoms()[atom_index + 1:], other_molecule.get_atoms()[atom_index + 1:]):

                # check to make sure that the atom1s have the same type
                if this_atom1.get_name() != other_atom1.get_name():
                    raise InconsistentValueError("self atom symbol", "other atom symbol", this_atom.get_name(), other_atom.get_name(), "symbols must be the same, make sure you are computing the rmsd of two molecules with the same atoms and fragments")

                # check to make sure that the atom2s have the same type
                if this_atom2.get_name() != other_atom2.get_name():
                    raise InconsistentValueError("self atom symbol", "other atom symbol", this_atom.get_name(), other_atom.get_name(), "symbols must be the same, make sure you are computing the rmsd of two molecules with the same atoms and fragments")

                # add these atom pairs' contribution to the squared distance difference
                squared_distance_difference += (this_atom1.distance(this_atom2) - other_atom1.distance(other_atom2)) ** 2

        # compute the rmsd of the sqrt of mean squared distance difference
        return math.sqrt(squared_distance_difference / self.get_num_atoms())

    def compare(self, other, cutoff_rmsd = 0.1):
        """
        Compares two molecules to see if they are similar to eachother bellow a cutoff rmsd

        Args:
            other - the molecule to compare this one to
            cutoff_rmsd - the rmsd level at which False will be returned, defailt is 0.1

        Returns:
            True if the rmsd between this molecule and the other is less than cutoff_rmsd, otherwise False
            Always returns False if the two molecules do not have the same fragments and atoms
        """

        try:
            return self.rmsd(other) < cutoff_rmsd
        except InconsistentValueError:
            return False

    def get_excluded_pairs(self, max_exclusion = 3):
        """
        Gets the excluded pairs of this molecule

        Args:
            None

        Returns:
            a tuple in the format (excluded_12, excluded_13, excluded_14, ..., excluded_1x) where each ecluded_1x is a list of lists of each fragment's excluded 1x pairs
        """

        excluded_pairs = [[] for i in range(max_exclusion)]

        for index, fragment in enumerate(self.get_fragments()):
            frag_excluded_pairs = fragment.get_excluded_pairs(max_exclusion)
            for exclusion_index in range(max_exclusion):
                excluded_pairs[exclusion_index].append(frag_excluded_pairs[exclusion_index])

        return excluded_pairs

    def to_xyz(self, fragments = None, cp = False):
        """
        Gets a string representation of the fragments in this molecule in the xyz file format

        Args:
            fragments   - list of fragment indicies to include in the string; optional, defualt is to include all fragments
            cp          - if True then fragments not specified in the fragments list will be included as ghost fragments

        Returns:
            String representation of the fragments in this molecule in the xyz format
        """

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

    def get_SHA1(self):
        """
        Generates the SHA1 hash of this molecule. Uses atoms, spin multiplicity and charge. Can be used to uniquely identify this molecule.

        Sorts fragments and atoms into standard order first, so the same molecule specified differently will have the same hash

        Args:
            None

        Returns:
            SHA1 hash of this molecule
        """

        hash_string = self.to_xyz() + "\n" + str(self.get_charge()) + "\n" + str(self.get_spin_multiplicity())
        return sha1(hash_string.encode()).hexdigest()

    def read_xyz(self, file, names, atoms_per_fragment, symmetry_per_fragment, charges, spin_multiplicities):
        """
        Reads a single xyz configuration from an open file into this molecule

        If the molecule already has some fragments, new fragments will be added in addition to any existing fragments

        Args:
            file  - the xyz file containing the atom coordinates
            names   - list of the names of each fragment
            atoms_per_fragment - list of the number of atoms in each fragment
            charges - list of the charges of each fragment
            spin_multiplicites - list of the spin multiplicities of each fragment

        Returns:
            This Molecule. Will be UNCHANGED if an empty file (or a file with all of its lines already read) is passed
        """

        first_line = file.readline()

        # check for end of file
        if first_line == "":
            raise StopIteration

        num_atoms = int(first_line)

        if not len(atoms_per_fragment) == len(symmetry_per_fragment):
            raise InconsistentValueError("atoms per fragment", "symmetry per fragment", atoms_per_fragment, symmetry_per_fragment, "lists must be same length")
        if not len(atoms_per_fragment) == len(charges):
            raise InconsistentValueError("atoms per fragment", "charges per fragment", atoms_per_fragment, charges, "lists must be same length")
        if not len(atoms_per_fragment) == len(spin_multiplicities):
            raise InconsistentValueError("atoms per fragment", "spin multiplicities per fragment", atoms_per_fragment, spin_multiplicities, "lists must be same length")
        if not len(atoms_per_fragment) == len(names):
            raise InconsistentValueError("atoms per fragment", "fragment names", atoms_per_fragment, names, "lists must be same length")

        if atoms_per_fragment != None and num_atoms != sum(atoms_per_fragment):
            raise InconsistentValueError("total atoms in xyz file", "fragments", num_atoms, atoms_per_fragment, "fragments list must sum to total atoms from input xyz file")

        # throw away the comment line
        file.readline()

        for name, atom_count, symmetry, charge, spin_multiplicity in zip(names, atoms_per_fragment, symmetry_per_fragment, charges, spin_multiplicities):
            # construct each fragment from its charge, spin multiplicity and its lines from the file
            self.add_fragment(Fragment(name, charge, spin_multiplicity).read_xyz(file, atom_count, symmetry))

        return self
