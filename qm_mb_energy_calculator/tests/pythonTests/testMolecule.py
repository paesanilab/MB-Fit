import sys
sys.path.insert(0, "../../src/")

import unittest
from molecule import Atom


"""
Test Cases for atoms
"""
class TestAtoms(unittest.TestCase):
    
    def test_name(self):
        atom = Atom("H", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_name(), "H")
        
        atom = Atom("Cl", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_name(), "Cl")
        
        atom = Atom("He", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_name(), "He")
        
        atom = Atom("Ar", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_name(), "Ar")

    def test_charge(self):
        atom = Atom("H", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_charge(), 0)
        
        atom = Atom("H", -3, 0, 0, 0, 0)
        self.assertEqual(atom.get_charge(), -3)
        
        atom = Atom("H", 7, 0, 0, 0, 0)
        self.assertEqual(atom.get_charge(), 7)
        
        atom = Atom("H", 1, 0, 0, 0, 0)
        self.assertEqual(atom.get_charge(), 1)

    def test_unpaired(self):
        atom = Atom("H", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_unpaired(), 0);
        
        atom = Atom("H", 0, 3, 0, 0, 0)
        self.assertEqual(atom.get_unpaired(), 3);
        
        atom = Atom("H", 0, 2, 0, 0, 0)
        self.assertEqual(atom.get_unpaired(), 2);
        
        atom = Atom("H", 0, 1, 0, 0, 0)
        self.assertEqual(atom.get_unpaired(), 1);

    def test_to_xyz(self):
        atom = Atom("H", 0, 0, 0, 0, 0)
        self.assertEqual(atom.to_xyz(), "H    0.00000000000000e+00   0.00000000000000e+00   0.00000000000000e+00");
        
        atom = Atom("Ar", 0, 0, 0.34234, -0.342, 1.2334)
        self.assertEqual(atom.to_xyz(), "Ar   3.42340000000000e-01  -3.42000000000000e-01   1.23340000000000e+00");
        
        atom = Atom("He", 0, 0, 12.34, 105.34, -0.00432)
        self.assertEqual(atom.to_xyz(), "He   1.23400000000000e+01   1.05340000000000e+02  -4.32000000000000e-03");
        
        atom = Atom("C", 0, 0, -2523452, 0.0003453, 34534)
        self.assertEqual(atom.to_xyz(), "C   -2.52345200000000e+06   3.45300000000000e-04   3.45340000000000e+04");

atomSuite = unittest.TestLoader().loadTestsFromTestCase(TestAtoms)

suite = unittest.TestSuite([atomSuite])
