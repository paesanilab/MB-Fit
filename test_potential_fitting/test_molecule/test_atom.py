import unittest

from potential_fitting.molecule import Atom

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
    def test_symmetry_class(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.get_symmetry_class(), "A")
        
        atom = Atom("Cl", "B", 0, 0, 0)
        self.assertEqual(atom.get_symmetry_class(), "B")
        
        atom = Atom("He", "C", 0, 0, 0)
        self.assertEqual(atom.get_symmetry_class(), "C")
        
        atom = Atom("Ar", "D", 0, 0, 0)
        self.assertEqual(atom.get_symmetry_class(), "D")

    """
    Tests the to_xyz() function of the Atom class
    """ 
    def test_to_xyz(self):
        atom = Atom("H", "A", 0, 0, 0)
        self.assertEqual(atom.to_xyz(), "H    0.00000000000000e+00   0.00000000000000e+00   0.00000000000000e+00");
        
        atom = Atom("Ar", "A", 0.34234, -0.342, 1.2334)
        self.assertEqual(atom.to_xyz(), "Ar   3.42340000000000e-01  -3.42000000000000e-01   1.23340000000000e+00");
        
        atom = Atom("He", "A", 12.34, 105.34, -0.00432)
        self.assertEqual(atom.to_xyz(), "He   1.23400000000000e+01   1.05340000000000e+02  -4.32000000000000e-03");
        
        atom = Atom("C", "A", -2523452, 0.0003453, 34534)
        self.assertEqual(atom.to_xyz(), "C   -2.52345200000000e+06   3.45300000000000e-04   3.45340000000000e+04");

        atom = Atom("Cl", "A", 0.00000000000068, 0.74576456847823583, 34534262462472457756745)
        self.assertEqual(atom.to_xyz(), "Cl   6.80000000000000e-13   7.45764568478236e-01   3.45342624624725e+22");

suite = unittest.TestLoader().loadTestsFromTestCase(TestAtom)
