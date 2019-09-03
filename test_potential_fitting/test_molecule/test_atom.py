import unittest, random, math

from potential_fitting.molecule import Atom
from potential_fitting.utils import Quaternion

"""
Test Cases for Atom class
"""
class TestAtom(unittest.TestCase):
   
    """
    Tests the get_name() function of the Atom class
    """ 
    def test_get_name(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.get_name(), "H")
        
        atom = Atom("Cl", "B", 0, 0, 0)
        self.assertEqual(atom.get_name(), "Cl")
        
        atom = Atom("He", "C", 0, 0, 0)
        self.assertEqual(atom.get_name(), "He")
        
        atom = Atom("Ar", "D", 0, 0, 0)
        self.assertEqual(atom.get_name(), "Ar")
 
    """
    Tests the get_symmetry_class() function of the Atom class
    """ 
    def test_get_symmetry_class(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.get_symmetry_class(), "A")
        
        atom = Atom("Cl", "B", 0, 0, 0)
        self.assertEqual(atom.get_symmetry_class(), "B")
        
        atom = Atom("He", "C", 0, 0, 0)
        self.assertEqual(atom.get_symmetry_class(), "C")
        
        atom = Atom("Ar", "D", 0, 0, 0)
        self.assertEqual(atom.get_symmetry_class(), "D")

    def test_set_symmetry_class(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.get_symmetry_class(), "A")

        atom.set_symmetry_class("B")
        self.assertEqual(atom.get_symmetry_class(), "B")

        atom.set_symmetry_class("C")
        self.assertEqual(atom.get_symmetry_class(), "C")

        atom.set_symmetry_class("A")
        self.assertEqual(atom.get_symmetry_class(), "A")

    def test_get_number(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.get_number(), 1)

        atom = Atom("Cl", "B", 0, 0, 0)
        self.assertEqual(atom.get_number(), 17)

        atom = Atom("He", "C", 0, 0, 0)
        self.assertEqual(atom.get_number(), 2)

        atom = Atom("Ar", "D", 0, 0, 0)
        self.assertEqual(atom.get_number(), 18)

    def test_get_mass(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.get_mass(), 1.008)

        atom = Atom("Cl", "B", 0, 0, 0)
        self.assertEqual(atom.get_mass(), 35.45)

        atom = Atom("He", "C", 0, 0, 0)
        self.assertEqual(atom.get_mass(), 4.0026)

        atom = Atom("Ar", "D", 0, 0, 0)
        self.assertEqual(atom.get_mass(), 39.948)

    def test_get_radius(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.get_radius(), 0.53)

        atom = Atom("Cl", "B", 0, 0, 0)
        self.assertEqual(atom.get_radius(), 0.79)

        atom = Atom("He", "C", 0, 0, 0)
        self.assertEqual(atom.get_radius(), 0.31)

        atom = Atom("Ar", "D", 0, 0, 0)
        self.assertEqual(atom.get_radius(), 0.71)

    def test_get_covalent_radius(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.get_covalent_radius(), 0.37)

        atom = Atom("Cl", "B", 0, 0, 0)
        self.assertEqual(atom.get_covalent_radius(), 0.99)

        atom = Atom("He", "C", 0, 0, 0)
        self.assertEqual(atom.get_covalent_radius(), 0.32)

        atom = Atom("Ar", "D", 0, 0, 0)
        self.assertEqual(atom.get_covalent_radius(), 0.97)

    def test_get_vdw_radius(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.get_vdw_radius(), 1.2)

        atom = Atom("Cl", "B", 0, 0, 0)
        self.assertEqual(atom.get_vdw_radius(), 1.75)

        atom = Atom("He", "C", 0, 0, 0)
        self.assertEqual(atom.get_vdw_radius(), 1.4)

        atom = Atom("Ar", "D", 0, 0, 0)
        self.assertEqual(atom.get_vdw_radius(), 1.88)

    def test_get_base_priority(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.get_base_priority(), atom.get_number())

        atom = Atom("Cl", "B", 0, 0, 0)
        self.assertEqual(atom.get_base_priority(), atom.get_number())

        atom = Atom("He", "C", 0, 0, 0)
        self.assertEqual(atom.get_base_priority(), atom.get_number())

        atom = Atom("Ar", "D", 0, 0, 0)
        self.assertEqual(atom.get_base_priority(), atom.get_number())

    def test_get_x_y_z(self):
        atom = Atom("H", "A", 0, 0, 0)

        self.assertEqual(atom.get_x(), 0)
        self.assertEqual(atom.get_y(), 0)
        self.assertEqual(atom.get_z(), 0)

        atom = Atom("Cl", "B", 0.76, 12.43, -436.234)

        self.assertEqual(atom.get_x(), 0.76)
        self.assertEqual(atom.get_y(), 12.43)
        self.assertEqual(atom.get_z(), -436.234)

    def test_set_x_y_z(self):
        atom = Atom("H", "A", 0, 0, 0)

        self.assertEqual(atom.get_x(), 0)
        self.assertEqual(atom.get_y(), 0)
        self.assertEqual(atom.get_z(), 0)

        atom.set_x(0.76)
        atom.set_y(12.43)
        atom.set_z(-436.234)

        self.assertEqual(atom.get_x(), 0.76)
        self.assertEqual(atom.get_y(), 12.43)
        self.assertEqual(atom.get_z(), -436.234)

    def test_set_xyz(self):
        atom = Atom("H", "A", 0, 0, 0)

        self.assertEqual(atom.get_x(), 0)
        self.assertEqual(atom.get_y(), 0)
        self.assertEqual(atom.get_z(), 0)

        atom.set_xyz(0.76, 12.43, -436.234)

        self.assertEqual(atom.get_x(), 0.76)
        self.assertEqual(atom.get_y(), 12.43)
        self.assertEqual(atom.get_z(), -436.234)

    def test_translate(self):
        atom = Atom("H", "A", 0, 0, 0)

        self.assertEqual(atom.get_x(), 0)
        self.assertEqual(atom.get_y(), 0)
        self.assertEqual(atom.get_z(), 0)

        atom.translate(0.76, 12.43, -436.234)

        self.assertEqual(atom.get_x(), 0.76)
        self.assertEqual(atom.get_y(), 12.43)
        self.assertEqual(atom.get_z(), -436.234)

        atom.translate(3, -3, 1.4)

        self.assertEqual(atom.get_x(), 3.76)
        self.assertEqual(atom.get_y(), 9.43)
        self.assertEqual(atom.get_z(), -434.834)

    def test_rotate(self):
        atom = Atom("H", "A", 0, 0, 0)

        for i in range(1000):

            ref_atom = Atom("H", "A", random.random() * 100, random.random() * 100, random.random() * 100)

            pre_dist = atom.distance(ref_atom)

            atom.rotate(Quaternion.get_random_rotation_quaternion(), origin_x=ref_atom.get_x(), origin_y=ref_atom.get_y(), origin_z=ref_atom.get_z())

            post_dist = atom.distance(ref_atom)

            self.assertAlmostEqual(pre_dist, post_dist)

    def test_distance(self):

        for i in range(1000):

            atom1 = Atom("H", "A", random.random() * 100, random.random() * 100, random.random() * 100)

            atom2 = Atom("Cl", "B", random.random() * 100, random.random() * 100, random.random() * 100)

            self.assertAlmostEqual(atom1.distance(atom2),
                                   math.sqrt((atom1.get_x() - atom2.get_x()) ** 2 +
                                             (atom1.get_y() - atom2.get_y()) ** 2 +
                                             (atom1.get_z() - atom2.get_z()) ** 2))

    def test_is_bonded(self):

        atom1 = Atom("H", "A", 0, 0, 0)

        atom2 = Atom("H", "A", 0.5, 0, 0)

        self.assertTrue(atom1.is_bonded(atom2, bond_sensitivity=1))

        atom1 = Atom("H", "A", 0, 0, 0)

        atom2 = Atom("H", "A", 1, 0, 0)

        self.assertFalse(atom1.is_bonded(atom2, bond_sensitivity=1))

        atom1 = Atom("H", "A", 0, 0, 0)

        atom2 = Atom("H", "A", 0.25, 0, 0)

        self.assertTrue(atom1.is_bonded(atom2, bond_sensitivity=0.5))

        atom1 = Atom("H", "A", 0, 0, 0)

        atom2 = Atom("H", "A", 0.5, 0, 0)

        self.assertFalse(atom1.is_bonded(atom2, bond_sensitivity=0.5))

    """
    Tests the to_xyz() function of the Atom class
    """

    def test_to_xyz(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.to_xyz(), "H    0.00000000000000e+00   0.00000000000000e+00   0.00000000000000e+00")

        atom = Atom("Ar", "A", 0.34234, -0.342, 1.2334)
        self.assertEqual(atom.to_xyz(), "Ar   3.42340000000000e-01  -3.42000000000000e-01   1.23340000000000e+00")

        atom = Atom("He", "A", 12.34, 105.34, -0.00432)
        self.assertEqual(atom.to_xyz(), "He   1.23400000000000e+01   1.05340000000000e+02  -4.32000000000000e-03")

        atom = Atom("C", "A", -2523452, 0.0003453, 34534)
        self.assertEqual(atom.to_xyz(), "C   -2.52345200000000e+06   3.45300000000000e-04   3.45340000000000e+04")

        atom = Atom("Cl", "A", 0.00000000000068, 0.74576456847823583, 34534262462472457756745)
        self.assertEqual(atom.to_xyz(), "Cl   6.80000000000000e-13   7.45764568478236e-01   3.45342624624725e+22")

    """
    Tests the to_ghost_xyz() function of the Atom class
    """
    def test_to_ghost_xyz(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.to_ghost_xyz(), "@H    0.00000000000000e+00   0.00000000000000e+00   0.00000000000000e+00")

        atom = Atom("Ar", "A", 0.34234, -0.342, 1.2334)
        self.assertEqual(atom.to_ghost_xyz(), "@Ar   3.42340000000000e-01  -3.42000000000000e-01   1.23340000000000e+00")

        atom = Atom("He", "A", 12.34, 105.34, -0.00432)
        self.assertEqual(atom.to_ghost_xyz(), "@He   1.23400000000000e+01   1.05340000000000e+02  -4.32000000000000e-03")

        atom = Atom("C", "A", -2523452, 0.0003453, 34534)
        self.assertEqual(atom.to_ghost_xyz(), "@C   -2.52345200000000e+06   3.45300000000000e-04   3.45340000000000e+04")

        atom = Atom("Cl", "A", 0.00000000000068, 0.74576456847823583, 34534262462472457756745)
        self.assertEqual(atom.to_ghost_xyz(), "@Cl   6.80000000000000e-13   7.45764568478236e-01   3.45342624624725e+22")

suite = unittest.TestLoader().loadTestsFromTestCase(TestAtom)
