import math, numpy as np
from potential_fitting.utils import constants
from potential_fitting.utils.math import test_difference_under_threshold

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
        self.set_xyz(x, y, z)

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

    def get_atomic_number(self):
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

    def get_vdw_radius(self):
        """
        Gets the vanderwalls radius of this atom

        Args:
            None

        Returns:
            The vanderwalls radius of this atom
        """
        return constants.symbol_to_vdw_radius(self.name)

    def get_base_priority(self):
        """
        Gets the base priority of this atom.
        This is equal to its atomic number.

        Args:
            None

        Returns:
            The priority of this atom.
        """

        return constants.symbol_to_number(self.name)


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

        if x == -0.0:
            x = 0.0
        self.x = x

    def set_y(self, y):
        """
        Sets the y position of this atom

        Args:
            y   - the new y position

        Returns:
            None
        """

        if y == -0.0:
            y = 0.0
        self.y = y

    def set_z(self, z):
        """
        Sets the z position of this atom

        Args:
            z   - the new z position

        Returns:
            None
        """

        if z == -0.0:
            z = 0.0
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

        self.set_x(x)
        self.set_y(y)
        self.set_z(z)

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

    def to_xyz(self, num_digits=14):
        """
        Gets the string representation of this atom in the xyz file format

        Args:
            num_digits - The number of digits after the decimal point to include when writing this atom's coordinates.
                    Default: 14 Maximum: 14

        Returns:
            String containing this atom's atomic symbol and coordinates in the xyz format
        """

        x = round(self.get_x(), num_digits)
        y = round(self.get_y(), num_digits)
        z = round(self.get_z(), num_digits)

        if x == -0.0:
            x = 0.0
        if y == -0.0:
            y = 0.0
        if z == -0.0:
            z = 0.0

        return "{0:2} {1:22.14e} {2:22.14e} {3:22.14e}".format(self.name, x, y, z)

    def is_bonded(self, atom, bond_sensitivity = 1.15):
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
    def to_ghost_xyz(self, num_digits=14):
        """
        Gets the string representation of this atom in the xyz file format as a ghost atom

        Args:
            num_digits - The number of digits after the decimal point to include when writing this atom's coordinates.
                    Default: 14 Maximum: 14

        Returns:
            String containing this atom's atomic symbol and coordinates in the xyz ghost atom format
        """

        return "@{}".format(self.to_xyz(num_digits=num_digits))

    def __eq__(self, other):
        return (self.get_name() == other.get_name() and self.get_symmetry_class() == other.get_symmetry_class()
                and test_difference_under_threshold(self.get_x(), other.get_x(), 0.00001) and test_difference_under_threshold(self.get_y(), other.get_y(), 0.00001) and test_difference_under_threshold(self.get_z(), other.get_z(), 0.00001))

    def __ne__(self, other):
        return not self == other
