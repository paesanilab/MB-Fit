import sys
sys.path.insert(0, "../../src/")

import unittest
from molecule import Atom
from molecule import Fragment
from molecule import Molecule


"""
Test Cases for Atom class
"""
class TestAtom(unittest.TestCase):
   
    """
    Tests the get_name() function of the Atom class
    """ 
    def test_get_name(self):
        atom = Atom("H", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_name(), "H")
        
        atom = Atom("Cl", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_name(), "Cl")
        
        atom = Atom("He", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_name(), "He")
        
        atom = Atom("Ar", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_name(), "Ar")

    """
    Tests the get_charge() function of the Atom class
    """ 
    def test_get_charge(self):
        atom = Atom("H", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_charge(), 0)
        
        atom = Atom("H", -3, 0, 0, 0, 0)
        self.assertEqual(atom.get_charge(), -3)
        
        atom = Atom("H", 7, 0, 0, 0, 0)
        self.assertEqual(atom.get_charge(), 7)
        
        atom = Atom("H", 1, 0, 0, 0, 0)
        self.assertEqual(atom.get_charge(), 1)

    """
    Tests the get_unpaired() function of the Atom class
    """ 
    def test_get_unpaired(self):
        atom = Atom("H", 0, 0, 0, 0, 0)
        self.assertEqual(atom.get_unpaired(), 0);
        
        atom = Atom("H", 0, 3, 0, 0, 0)
        self.assertEqual(atom.get_unpaired(), 3);
        
        atom = Atom("H", 0, 2, 0, 0, 0)
        self.assertEqual(atom.get_unpaired(), 2);
        
        atom = Atom("H", 0, 1, 0, 0, 0)
        self.assertEqual(atom.get_unpaired(), 1);

    """
    Tests the to_xyz() function of the Atom class
    """ 
    def test_to_xyz(self):
        atom = Atom("H", 0, 0, 0, 0, 0)
        self.assertEqual(atom.to_xyz(), "H    0.00000000000000e+00   0.00000000000000e+00   0.00000000000000e+00");
        
        atom = Atom("Ar", 0, 0, 0.34234, -0.342, 1.2334)
        self.assertEqual(atom.to_xyz(), "Ar   3.42340000000000e-01  -3.42000000000000e-01   1.23340000000000e+00");
        
        atom = Atom("He", 0, 0, 12.34, 105.34, -0.00432)
        self.assertEqual(atom.to_xyz(), "He   1.23400000000000e+01   1.05340000000000e+02  -4.32000000000000e-03");
        
        atom = Atom("C", 0, 0, -2523452, 0.0003453, 34534)
        self.assertEqual(atom.to_xyz(), "C   -2.52345200000000e+06   3.45300000000000e-04   3.45340000000000e+04");

"""
Test cases for Fragment class
"""
class TestFragment(unittest.TestCase):
    

    """
    Tests the add_atom() and get_atom() functions of the Fragment class
    """
    def test_add_atom_and_get_atom(self):
        fragment = Fragment()

        atom0 = Atom("H", 0, 0, 0, 0, 0)

        fragment.add_atom(atom0)

        self.assertEqual(fragment.get_atoms()[0], atom0)

        atom1 = Atom("Cl", 1, 2, 400, 32, 23)    
    
        fragment.add_atom(atom1)
        
        self.assertEqual(fragment.get_atoms()[0], atom0)
        self.assertEqual(fragment.get_atoms()[1], atom1)

    """
    Tests the get_charge() function of the Fragment class
    """
    def test_get_charge(self):
        fragment = Fragment()
        
        # get_charge() should be 0 before any Atoms added to Fragment
        self.assertEqual(fragment.get_charge(), 0)

        fragment.add_atom(Atom("H", 0, 0, 0, 0, 0))
        
        # get_charge() should still be 0 after Atom with 0 charge added to Fragment
        self.assertEqual(fragment.get_charge(), 0)

        fragment.add_atom(Atom("Cl", -1, 0, 0, 0, 0))

        # get_charge() should be -1 after clorine ion added to Fragment
        self.assertEqual(fragment.get_charge(), -1)

        fragment.add_atom(Atom("F", -1, 0, 0, 0, 0))

        # get_charge() should be -2 after flourine ion added to Fragment
        self.assertEqual(fragment.get_charge(), -2)

        fragment.add_atom(Atom("K", 1, 0, 0, 0, 0))

        # get_charge() should be -1 after Potassium ion added to Fragment
        self.assertEqual(fragment.get_charge(), -1)

    """
    Tests the get_unpaired() function of the Fragment class
    """
    def test_get_unpaired(self):
        fragment = Fragment();

        # get_unpaired() should be 0 before any Atoms added to Fragment
        self.assertEqual(fragment.get_unpaired(), 0)

        fragment.add_atom(Atom("H", 0, 1, 0, 0, 0))

        # get_unpaired() should be 1 after Hydrogen added to fragment
        self.assertEqual(fragment.get_unpaired(), 1)

        fragment.add_atom(Atom("Tc", 0, 5, 0, 0, 0))

        # get_unpaired() should be 6 after Technetium added to fragment
        self.assertEqual(fragment.get_unpaired(), 6)



atomSuite = unittest.TestLoader().loadTestsFromTestCase(TestAtom)
fragmentSuite = unittest.TestLoader().loadTestsFromTestCase(TestFragment)

suite = unittest.TestSuite([atomSuite, fragmentSuite])
