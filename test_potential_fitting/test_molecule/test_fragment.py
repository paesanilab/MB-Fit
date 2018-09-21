import unittest

from potential_fitting.molecule import Atom
from potential_fitting.molecule import Fragment

"""
Test cases for Fragment class
"""
class TestFragment(unittest.TestCase):
    

    """
    Tests the add_atom() and get_atoms() functions of the Fragment class
    """
    def test_add_atom_and_get_atoms(self):
        fragment = Fragment("HCl", 0, 1)

        # get_atom() should return list of length 0 before any atoms added to fragment
        self.assertEqual(len(fragment.get_atoms()), 0)

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment.add_atom(atom0)

        self.assertEqual(fragment.get_atoms()[0], atom0)
        # get_atom() should return list of length 1 after 1 atom added to fragment
        self.assertEqual(len(fragment.get_atoms()), 1)

        atom1 = Atom("Cl", "B", 400, 32, 23)    
    
        fragment.add_atom(atom1)
        
        self.assertEqual(fragment.get_atoms()[0], atom0)
        self.assertEqual(fragment.get_atoms()[1], atom1)
        # get_atom() should return list of length 2 after 2 atoms added to fragment
        self.assertEqual(len(fragment.get_atoms()), 2)

    """
    Tests the get_charge() function of the Fragment class
    """
    def test_get_charge(self):
        fragment = Fragment("Name", 0, 2)
        
        self.assertEqual(fragment.get_charge(), 0)        

        fragment = Fragment("Name", -3, 1)
        
        self.assertEqual(fragment.get_charge(), -3)        
        
        fragment = Fragment("Name", 12, 5)
        
        self.assertEqual(fragment.get_charge(), 12)        


    """
    Tests the get_spin_multiplicity() function of the Fragment class
    """
    def test_get_spin_multiplicity(self):
        fragment = Fragment("Name", 0, 2)
        
        self.assertEqual(fragment.get_spin_multiplicity(), 2)        

        fragment = Fragment("Name", -3, 1)
        
        self.assertEqual(fragment.get_spin_multiplicity(), 1)        
        
        fragment = Fragment("Name", 12, 5)
        
        self.assertEqual(fragment.get_spin_multiplicity(), 5)        


    """
    Tests the get_num_atoms() function of the Fragment class
    """
    def test_get_num_atoms(self):
        fragment = Fragment("HAlHe", 3, 2);

        # get_num_atoms() should return 0 before any Atoms added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 0)
        
        fragment.add_atom(Atom("H", "A", 0, 0, 0))

        # get_num_atoms() should return 1 after single atom added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 1)
        
        fragment.add_atom(Atom("Al", "B", 0, 0, 0))

        # get_num_atoms() should return 2 after second atom added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 2)
        
        fragment.add_atom(Atom("He", "C", 0, 0, 0))

        # get_num_atoms() should return 3 after third atom added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 3)

    """
    Test the to_xyz() function of the Fragment class
    """
    def test_to_xyz(self):
        fragment = Fragment("HClXe", -2, 2);

        # to_xyz() should return empty string when no atoms are added to fragment
        self.assertEqual(fragment.to_xyz(), "");
        
        atom0 = Atom("H", "A", 0, 0, 0)

        fragment.add_atom(atom0)

        # to_xyz() should return string of first atom after only 1 atom added
        self.assertEqual(fragment.to_xyz(), atom0.to_xyz() + "\n")
        
        atom1 = Atom("Cl", "B", 5, 7, -3)

        fragment.add_atom(atom1)

        # to_xyz() should return string of 2 atoms after 2nd atom added
        self.assertEqual(fragment.to_xyz(), atom0.to_xyz() + "\n" + atom1.to_xyz() + "\n")
    
        atom2 = Atom("Xe", "C", 10.234235, -0.00000234, 2.353523)

        fragment.add_atom(atom2)

        # to_xyz() should return string of 3 atoms after only 3rd atom added
        self.assertEqual(fragment.to_xyz(), atom0.to_xyz() + "\n" + atom1.to_xyz() + "\n" + atom2.to_xyz() + "\n")

suite = unittest.TestLoader().loadTestsFromTestCase(TestFragment)
