import numpy, math

from hashlib import sha1

from potential_fitting.exceptions import XYZFormatError, InvalidValueError, InconsistentValueError
from .fragment import Fragment

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
        try:
            symmetry = self.get_fragments()[0].get_symmetry()
        except IndexError:

            # if there are no fragments, symmetry is empty string
            return ""

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

        fragment.index = self.get_num_fragments()

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

        return sorted(self.fragments, key=lambda fragment: fragment.get_index())

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

    def get_symbols(self):
        """
        Gets the atomic symbols of the atoms in this molecule as a list

        Args:
            None

        Returns:
            list of the atomic symbols of the atoms in this molecule
        """

        return [atom.get_name() for atom in self.get_atoms()]

    def get_coordinates(self):
        """
        Gets the positions of the atoms in this molecule as a list of 3-tuples

        Args:
            None

        Returns:
            list of the positions of the atoms in this moleule
        """

        return [(atom.get_x(), atom.get_y(), atom.get_z()) for atom in self.get_atoms()]

    def read_fragment_from_xyz(self, string, name, charge, spin_multiplicity, symmetry):
        """
        Reads a single fragment from the given string and adds it to this molecule

        Args:
            string - the string in the xyz file format, should include just the atom lines, no total atom line or comment line
            name - the name of this new fragment
            charge - the charge of this new fragment
            spin_multiplicity - the spin multiplicity of the new fragment
            symmetry - the symmetry of the new fragment, specified as a string in form A1B2

        Returns:
            self
        """

        self.add_fragment(Fragment(name, charge, spin_multiplicity).read_xyz(string, symmetry))

        return self

    def read_xyz(self, string, atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_multiplicity_per_fragment, symmetry_per_fragment):
        """
        Reads fragments from an xyz string into this Molecule

        Args:
            string - the xyz format string
            atoms_per_fragment - list containing the number of atoms in each fragment
            name_per_fragment - list containing the names of each fragment
            charge_per_fragment - list containing the charges of each fragment
            spin_multiplicity_per_fragment - list containing the spin multiplicities of each fragment
            symmetry_per_fragment - list containing the symmetries of each fragment, in format A1B2

        Returns:
            self
        """

        # Error checking to make sure all lists passed in are the same length
        if not len(atoms_per_fragment) == len(symmetry_per_fragment):
            raise InconsistentValueError("atoms per fragment", "symmetry per fragment", atoms_per_fragment, symmetry_per_fragment, "lists must be same length")
        if not len(atoms_per_fragment) == len(charge_per_fragment):
            raise InconsistentValueError("atoms per fragment", "charges per fragment", atoms_per_fragment, charge_per_fragment, "lists must be same length")
        if not len(atoms_per_fragment) == len(spin_multiplicity_per_fragment):
            raise InconsistentValueError("atoms per fragment", "spin multiplicities per fragment", atoms_per_fragment, spin_multiplicity_per_fragment, "lists must be same length")
        if not len(atoms_per_fragment) == len(name_per_fragment):
            raise InconsistentValueError("atoms per fragment", "fragment names", atoms_per_fragment, names, "lists must be same length")

        # break the input string apart along \n characters
        lines = string.splitlines()

        # read the total number of atoms from the first line of the xyz
        try:
            atom_total = int(lines[0])
        except ValueError:
            raise XYZFormatError("atom count line '{}' cannot be parsed into an integer".format(lines[0]), "line should contain a single integer")

        # make sure that the total number of atoms indicated by the xyz file matches the number of atoms indicated per fragment
        if atom_total != sum(atoms_per_fragment):
            raise InconsistentValueError("total atoms in xyz string", "fragments", atom_total, atoms_per_fragment, "fragments list must sum to total atoms from input xyz string")

        # remove the atom total and comment lines from the lines list
        lines = lines[2:]

        # make sure that there are a number of lines equal to the total number of atoms
        if len(lines) != atom_total:
            raise InconsistentValueError("total atoms in xyz string", "atom lines in xyz string", atom_total, len(lines), "number of total atoms indicated in xyz string should match number of atom lines")

        # loop over each item in the lists, each iteration containing the information to assemble one fragment
        for num_atoms, name, charge, spin, symmetry in zip(atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_multiplicity_per_fragment, symmetry_per_fragment):

            # read this fragment into this Molecule
            self.read_fragment_from_xyz("\n".join(lines[:num_atoms]), name, charge, spin, symmetry)

            # remove a number of lines from the lines list equal to the number used in the Fragment that was just read
            lines = lines[num_atoms:]

        return self

    def read_xyz_file(self, file, atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_multiplicity_per_fragment, symmetry_per_fragment):
        """
        Reads fragments from an xyz file into this Molecule

        Will attempt to read lines from the given file handle, raising a StopIteration exception if called on an empty file and a
        

        Args:
            file - the file to read from
            atoms_per_fragment - list containing the number of atoms in each fragment
            name_per_fragment - list containing the names of each fragment
            charge_per_fragment - list containing the charges of each fragment
            spin_multiplicity_per_fragment - list containing the spin multiplicities of each fragment
            symmetry_per_fragment - list containing the symmetries of each fragment, in format A1B2

        Returns:
            self
        """
        
        # build the xyz string
        string = ""

        # read lines from the file equal to the number needed for one molecule
        for line_count in range(2 + sum(atoms_per_fragment)):

            line = file.readline()

            # if the line is an empty string, then we have reached end of file mid parse
            if line == "":
                if line_count == 0:
                    raise StopIteration # if the first line is empty, raise StopIteration to indicate that this file is out of molecules to parse
                raise XYZFormatError("ran out of lines to read from xyz file {} in the middle of a molecule".format(file.name), "make sure the last molecule in the file has a comment line and a number of atoms equal to the amount indicated in the atom count line.")

            string += line
        
        self.read_xyz(string, atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_multiplicity_per_fragment, symmetry_per_fragment)

        return self

    def read_xyz_path(self, path, atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_multiplicity_per_fragment, symmetry_per_fragment):
        """
        Reads fragments from an xyz file indicated by a filepath into this Molecule

        Will attempt to read lines from the file at the given file path, raising an exception if it runs out of lines mid-parse

        Args:
            path - the path to the file to read from
            atoms_per_fragment - list containing the number of atoms in each fragment
            name_per_fragment - list containing the names of each fragment
            charge_per_fragment - list containing the charges of each fragment
            spin_multiplicity_per_fragment - list containing the spin multiplicities of each fragment
            symmetry_per_fragment - list containing the symmetries of each fragment, in format A1B2

        Returns:
            self
        """

        with open(path, "r") as file:

            try:
                self.read_xyz_file(file, atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_multiplicity_per_fragment, symmetry_per_fragment)

            # if the call to read_xyz_file() raises a StopIteration, it means the file was empty
            except StopIteration:
                raise XYZFormatError("xyz file {} file is empty".format(file.name), "make sure the xyz file has at least 1 molecule in it")

        return self

    def read_xyz_direct(self, string, settings = None):
        """
        Reads fragments from a string into this Molecule

        Will infer a single fragment with charge 0, spin 1, and no symmetry if settings is None

        Args:
            string - the string to read from
            settings - settings file containing information about the molecule

        Returns:
            self
        """

        # if settings is None, then infer default values for molecule attributes
        if settings is None:
            name_per_fragment = ["noname"]
            charge_per_fragment = [0]
            spin_multiplicity_per_fragment = [1]

            total_atoms = int(string.splitlines()[0])

            atoms_per_fragment = [total_atoms]        

            symmetry = ""

            try:
                symmetry_class = ord(max(self.get_symmetry().strip("_"))) + 1
                print("SYM START:", symmetry) 
            except ValueError:
                symmetry_class = 65

            # loop over each atom assigning it a unique symmetry class
            for atom_index in range(total_atoms):

                symmetry += "{}1".format(chr(symmetry_class))

                symmetry_class += 1


            symmetry_per_fragment = [symmetry]
            
        # if settings is defined, read values from xyz file
        else:
            atoms_per_fragment = [int(count) for count in settings.get("molecule", "fragments").split(",")]
            name_per_fragment = [int(name) for name in settings.get("molecule", "name").split(",")]
            charge_per_fragment = [int(charge) for charge in settings.get("molecule", "charges").split(",")]
            spin_multiplicity_per_fragment = [int(spin) for spin in settings.get("molecule", "spin").split(",")]
            symmetry_per_fragment = [int(symmetry) for symmetry in settings.get("molecule", "symmetry").split(",")]


        self.read_xyz(string, atoms_per_fragment, name_per_fragment, charge_per_fragment, spin_multiplicity_per_fragment, symmetry_per_fragment)

        return self

    def read_xyz_file_direct(self, file, settings = None): 
        """
        Reads fragments from a file into this Molecule

        Will infer a single fragment with charge 0, spin 1, and no symmetry if settings is None

        Args:
            string - the string to read from
            settings - settings file containing information about the molecule

        Returns:
            self
        """
        if settings is None:
            position = file.tell()
            atoms_per_fragment = [int(file.readline())]
            file.seek(position)
        else:
            atoms_per_fragment = [int(count) for count in settings.get("molecule", "fragments").split(",")]
        # build the xyz string
        string = ""

        # read lines from the file equal to the number needed for one molecule
        for line_count in range(2 + sum(atoms_per_fragment)):

            line = file.readline()

            # if the line is an empty string, then we have reached end of file mid parse
            if line == "":
                if line_count == 0:
                    raise StopIteration # if the first line is empty, raise StopIteration to indicate that this file is out of molecules to parse
                raise XYZFormatError("ran out of lines to read from xyz file {} in the middle of a molecule".format(file.name), "make sure the last molecule in the file has a comment line and a number of atoms equal to the amount indicated in the atom count line.")

            string += line
        
        self.read_xyz_direct(string, settings)

        return self

    def read_xyz_path_direct(self, path, settings = None):
        """
        Reads fragments from an xyz_file indicated by a path into this Molecule

        Will infer a single fragment with charge 0, spin 1, and no symmetry if settings is None

        Args:
            string - the string to read from
            settings - settings file containing information about the molecule

        Returns:
            self
        """
        
        with open(path, "r") as file:

            try:
                self.read_xyz_file_direct(file, settings)

            # if the call to read_xyz_file() raises a StopIteration, it means the file was empty
            except StopIteration:
                raise XYZFormatError("xyz file {} file is empty".format(file.name), "make sure the xyz file has at least 1 molecule in it")

        return self

    def read_psi4_string(self, string):
        """
        Reads the string outputted by a call to psi4.molecule.save_string_xyz() into this molecule as a fragment

        the fragments added will not have name or symmetry saved correctly, because this information is not available
        from the output of psi4.molecule.save_string_xyz(). As a result certain operations will not work on this molecule, for example
        do not add this molecule to a database or attempt to generate its polynomial input format in style A1B2

        Should never be called on a molecule that already has fragments, as this will have unexpected results

        Args:
            output of psi4.molecule.save_string_xyz()

        Returns:
            self
        """

        # divide the string along \n characters
        lines = string.splitlines()

        # read charge and spin from first line of input string, casting each to an int
        try:
            charge, spin_multiplicity = [int(value) for value in lines[0].split()]
        except ValueError:
            raise XYZFormatError(lines[0], "line format should be 'charge spin_multiplicity', make sure you are passing in the output of psi4.molecule.save_string_xyz()")

        # calculate total atoms in this molecule
        total_atoms = len(lines) - 1

        # these fields do not matter
        name = "unnamed"

        # used to build the symmetry string for the fragment
        symmetry = ""

        # keeps track of which symmetry_class to use for the next atom
        symmetry_class = 65

        # loop over each atom assigning it a unique symmetry class
        for atom_index in range(total_atoms):

            symmetry += "{}1".format(chr(symmetry_class))

            symmetry_class += 1

        self.read_fragment_from_xyz("\n".join(lines[1:]), name, charge, spin_multiplicity, symmetry)

        return self
