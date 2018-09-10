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

    def __init__(self, name, x, y, z):
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

    def rotate(self, x_radians, y_radians, z_radians, x_origin = 0, y_origin = 0, z_origin = 0):
        """
        Rotates this atom around a point

        Args:
            x_radians - the number of radians to rotate around the x axis
            y_radians - the number of radians to rotate around the y axis
            z_radians - the number of radians to rotate around the z axis
            x_origin - x position of point to rotate around, default is 0
            y_origin - y position of point to rotate around, default is 0
            z_origin - z position of point to rotate around, default is 0

        Returns:
            None
        """

        # first construct the matrix of rotation
        rotation_x_matrix = numpy.matrix([
            [1,                         0,                      0                       ],
            [0,                         math.cos(x_radians),    - math.sin(x_radians)   ], 
            [0,                         math.sin(x_radians),    math.cos(x_radians)     ]
                                        ])

        rotation_y_matrix = numpy.matrix([
            [math.cos(y_radians),       0,                      math.sin(y_radians)     ],
            [0,                         1,                      0                       ], 
            [- math.sin(y_radians),     0,                      math.cos(y_radians)     ]
                                        ])

        rotation_z_matrix = numpy.matrix([
            [math.cos(z_radians),       - math.sin(z_radians),  0                       ],
            [math.sin(z_radians),       math.cos(z_radians),    0                       ], 
            [0,                         0,                      1                       ]
                                        ])

        rotation_matrix = rotation_x_matrix * rotation_y_matrix * rotation_z_matrix

        # get the new xyz values after multiplying by rotation matrix, moving coordinates to transform around the origin
        x, y, z = (numpy.matrix([self.x - x_origin, self.y - y_origin, self.z - z_origin]) * numpy.matrix(rotation_matrix)).getA1()

        # update the xyz position of this atom, adding back the origin coordinates
        self.set_xyz(x + x_origin, y + y_origin, z + z_origin)
        
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

        self.atoms.append(atom)

    def get_atoms(self):
        """
        Gets a list of the atoms in this fragment, sorted in standard order (alphabetically by xyz representation

        Args:
            None

        Returns:
            A list of the atoms in this fragment, sorted in standard order
        """

        return sorted(self.atoms, key=lambda atom: atom.to_xyz())

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

    def rotate(self, x_radians, y_radians, z_radians, x_origin = 0, y_origin = 0, z_origin = 0):
        """
        Rotates this fragment around a point

        Args:
            x_radians - the number of radians to rotate around the x axis
            y_radians - the number of radians to rotate around the y axis
            z_radians - the number of radians to rotate around the z axis
            x_origin - x position of point to rotate around, default is 0
            y_origin - y position of point to rotate around, default is 0
            z_origin - z position of point to rotate around, default is 0

        Returns:
            None
        """

        for atom in self.get_atoms():
            atom.rotate(x_radians, y_radians, z_radians, x_origin, y_origin, z_origin)

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
    
    def read_xyz(self, file, num_atoms):
        """
        Reads a number of lines from an open file into this fragment

        Args:
            file  - the file to read lines from
            num_atoms - the number of lines (and atoms) to read from the xyz file

        Returns:
            This Fragment
        """

        # loop once for each atom expected to be read
        for i in range(num_atoms):

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

            self.add_atom(Atom(symbol, float(x), float(y), float(z)))

        return self

    def read_xyz_string(self, string):
        """
        Reads the given string's atoms into this fragment

        Args:
            string - the xyz file string

        Returns:
            None
        """

        for line in string.splitlines():
            try:
                symbol, x, y, z = line.split()
            except ValueError:
                raise XYZFormatError(line, "ATOMIC_SYMBOL X Y Z") from None

            self.add_atom(Atom(symbol, float(x), float(y), float(z)))

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

    def add_fragment(self, fragment):
        """
        Adds a fragment to this molecule

        Args:
            fragment    - the fragment to add

        Returns:
            None
        """

        self.fragments.append(fragment)

    def get_fragments(self):
        """
        Gets a list of the fragments in this molecule in standard order

        Args:
            None

        Returns:
            List of fragments in this molecule in standard order
        """

        return sorted(self.fragments, key=lambda fragment: fragment.get_name() + " " + fragment.to_xyz())

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

    def rotate(self, x_radians, y_radians, z_radians, x_origin = 0, y_origin = 0, z_origin = 0):
        """
        Rotates this molecule around a point

        Args:
            x_radians - the number of radians to rotate around the x axis
            y_radians - the number of radians to rotate around the y axis
            z_radians - the number of radians to rotate around the z axis
            x_origin - x position of point to rotate around, default is 0
            y_origin - y position of point to rotate around, default is 0
            z_origin - z position of point to rotate around, default is 0

        Returns:
            None
        """

        for fragment in self.get_fragments():
            fragment.rotate(x_radians, y_radians, z_radians, x_origin, y_origin, z_origin)

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
            atom.set_xyz(x, y, z)

    def rmsd(self, other):
        """
        Computes the RMSD distance between the atoms in two molecules

        molecules must have the same fragments and atoms or an InconsistentValueError will be raised

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

    def read_xyz(self, file, names, atoms_per_fragment, charges, spin_multiplicities):
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

        if not len(atoms_per_fragment) == len(charges):
            raise InconsistentValueError("atoms per fragment", "charges per fragment", atoms_per_fragment, charges, "lists must be same length")
        if not len(atoms_per_fragment) == len(spin_multiplicities):
            raise InconsistentValueError("atoms per fragment", "spin multiplicities per fragment", atoms_per_fragment, spin_multiplicities, "lists must be same length")
        if not len(atoms_per_fragment) == len(names):
            raise InconsistentValueError("atoms per fragment", "fragment names", atoms_per_fragment, names, "lists must be same length")

        if num_atoms != sum(atoms_per_fragment):
            raise InconsistentValueError("total atoms in xyz file", "fragments", num_atoms, atoms_per_fragment, "fragments list must sum to total atoms from input xyz file")

        # throw away the comment line
        file.readline()

        for name, atom_count, charge, spin_multiplicity in zip(names, atoms_per_fragment, charges, spin_multiplicities):
            # construct each fragment from its charge, spin multiplicity and its line's from the string
            self.add_fragment(Fragment(name, charge, spin_multiplicity).read_xyz(file, atom_count))

        return self

    def read_psi4_string(self, string, name):
        """
        Reads an xyz file format with 
        """
        lines = string.splitlines();

        charge, spin = lines[0].split()

        self.add_fragment(Fragment(name, float(charge), int(spin)).read_xyz_string("\n".join(lines[1:])))

        return self
