import unittest
import os

from potential_fitting.molecule import Atom
from potential_fitting.molecule import Fragment
from potential_fitting.molecule import parse_training_set_file
from potential_fitting.utils import SettingsReader
from potential_fitting.exceptions import InvalidValueError, InconsistentValueError

"""
Test cases for Fragment class
"""
class TestFragment(unittest.TestCase):

    def test_bad_SMILE(self):

        # tests for when a fragment is given a bad SMILE string

        # SMILE string that doesn't close all its bonds.
        with self.assertRaises(InvalidValueError):
            Fragment([Atom("O", "A", 0, 0, 0)], "frag", 0, 1, "O1")
        with self.assertRaises(InvalidValueError):
            Fragment([Atom("H", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("C", "C", 2, 2, 2)], "frag", 0, 1, "O12HC2")

        # SMILE string that indicates an atom is bonded to itself
        with self.assertRaises(InvalidValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("H", "B", 2, 2, 2)], "frag", 0, 1, "O11(H)H")

    def test_inconsistent_atoms_and_SMILE(self):

        # tests for when atoms in the list don't match the atoms in the SMILE

        # too few atoms in list
        with self.assertRaises(InconsistentValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1)], "frag", 0, 1, "O(H)H")

        # too many atoms in list
        with self.assertRaises(InconsistentValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("H", "B", 2, 2, 2), Atom("H", "B", 3, 3, 3)],
                "frag", 0, 1, "O(H)H")

        # inconsistent atom names
        with self.assertRaises(InconsistentValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("O", "A", 1, 1, 1), Atom("H", "B", 2, 2, 2)], "frag", 0, 1, "O(H)H")

    def test_negative_spin(self):

        # test for when a fragment is given a spin multiplicity < 1

        with self.assertRaises(InvalidValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("H", "B", 2, 2, 2)], "frag", 0, 0, "O(H)H")

        with self.assertRaises(InvalidValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("H", "B", 2, 2, 2)], "frag", 0, -1, "O(H)H")

    def test_get_name(self):
        fragment = Fragment([], "HCl", 0, 1, "")
        self.assertEqual(fragment.get_name(), "HCl")

        fragment = Fragment([], "somename!", 0, 1, "")
        self.assertEqual(fragment.get_name(), "somename!")

    def test_set_name(self):
        fragment = Fragment([], "HCl", 0, 1, "")
        self.assertEqual(fragment.get_name(), "HCl")

        fragment.set_name("somename!")
        self.assertEqual(fragment.get_name(), "somename!")

        fragment.set_name("another name")
        self.assertEqual(fragment.get_name(), "another name")

        fragment.set_name("HCl")
        self.assertEqual(fragment.get_name(), "HCl")

    """
    Tests get_atoms() functions of the Fragment class
    """
    def test_get_atoms(self):
        fragment = Fragment([], "HCl", 0, 1, "")

        # get_atom() should return list of length 0 before any atoms added to fragment
        self.assertEqual(len(fragment.get_atoms()), 0)

        atom0 = Atom("H", "A", 0, 0, 0)
        fragment = Fragment([atom0], "HCl", 0, 1, "H")

        self.assertEqual(fragment.get_atoms()[0], atom0)
        # get_atom() should return list of length 1 after 1 atom added to fragment
        self.assertEqual(len(fragment.get_atoms()), 1)

        atom1 = Atom("Cl", "B", 400, 32, 23)    

        fragment = Fragment([atom0, atom1], "HCl", 0, 1, "H[Cl]")
        
        self.assertEqual(fragment.get_atoms()[0], atom0)
        self.assertEqual(fragment.get_atoms()[1], atom1)

        # get_atom() should return list of length 2 after 2 atoms added to fragment
        self.assertEqual(len(fragment.get_atoms()), 2)

    """
    Tests the get_charge() function of the Fragment class
    """
    def test_get_charge(self):
        fragment = Fragment([], "Name", 0, 2, "")
        
        self.assertEqual(fragment.get_charge(), 0)        

        fragment = Fragment([], "Name", -3, 1, "")
        
        self.assertEqual(fragment.get_charge(), -3)        
        
        fragment = Fragment([], "Name", 12, 5, "")
        
        self.assertEqual(fragment.get_charge(), 12)        


    """
    Tests the get_spin_multiplicity() function of the Fragment class
    """
    def test_get_spin_multiplicity(self):
        fragment = Fragment([], "Name", 0, 2, "")
        
        self.assertEqual(fragment.get_spin_multiplicity(), 2)        

        fragment = Fragment([], "Name", -3, 1, "")
        
        self.assertEqual(fragment.get_spin_multiplicity(), 1)        
        
        fragment = Fragment([], "Name", 12, 5, "")
        
        self.assertEqual(fragment.get_spin_multiplicity(), 5)        


    """
    Tests the get_num_atoms() function of the Fragment class
    """
    def test_get_num_atoms(self):
        fragment = Fragment([], "HAlHe", 3, 2, "")

        # get_num_atoms() should return 0 before any Atoms added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 0)

        fragment = Fragment([Atom("H", "A", 0, 0, 0)], "HAlHe", 3, 2, "H")

        # get_num_atoms() should return 1 after single atom added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 1)

        fragment = Fragment([Atom("H", "A", 0, 0, 0), Atom("Al", "B", 0, 0, 0)], "HAlHe", 3, 2, "H[Al]")

        # get_num_atoms() should return 2 after second atom added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 2)

        fragment = Fragment([Atom("H", "A", 0, 0, 0), Atom("Al", "B", 0, 0, 0), Atom("He", "C", 0, 0, 0)], "HAlHe", 3, 2, "H[Al][He]")

        # get_num_atoms() should return 3 after third atom added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 3)

    """
    Test the to_xyz() function of the Fragment class
    """
    def test_to_xyz(self):
        fragment = Fragment([], "HClXe", -2, 2, "")

        # to_xyz() should return empty string when no atoms are added to fragment
        self.assertEqual(fragment.to_xyz(), "")

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "HClXe", -2, 2, "H")

        # to_xyz() should return string of first atom after only 1 atom added
        self.assertEqual(fragment.to_xyz(), atom0.to_xyz() + "\n")
        
        atom1 = Atom("Cl", "B", 5, 7, -3)

        fragment = Fragment([atom0, atom1], "HClXe", -2, 2, "H[Cl]")

        # to_xyz() should return string of 2 atoms after 2nd atom added
        self.assertEqual(fragment.to_xyz(), atom0.to_xyz() + "\n" + atom1.to_xyz() + "\n")
    
        atom2 = Atom("Xe", "C", 10.234235, -0.00000234, 2.353523)

        fragment = Fragment([atom0, atom1, atom2], "HClXe", -2, 2, "H[Cl][Xe]")

        # to_xyz() should return string of 3 atoms after only 3rd atom added
        self.assertEqual(fragment.to_xyz(), atom0.to_xyz() + "\n" + atom1.to_xyz() + "\n" + atom2.to_xyz() + "\n")

    """
    Test the get_SMILE() function of the Fragment class
    """
    def test_get_SMILE(self):
        # A4B2C2D4E4
        fragment1 = list(parse_training_set_file(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "bdc.xyz"), SettingsReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "bdc.ini"))))[0].get_fragments()[0]
        self.assertEqual(fragment1.get_SMILE(), "[C]%1%2[C]%3%4.[C]%5%6[C]%7%8.[C]%3%5%9.[C]%1%7%10.[H]%2.[H]%4.[H]%6.[H]%8.[C]%9%11%12.[C]%10%13%14.[O]%11.[O]%12.[O]%13.[O]%14")

        # A4B2C2D4E4
        fragment2 = list(parse_training_set_file(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "bdc2.xyz"), SettingsReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "bdc2.ini"))))[0].get_fragments()[0]
        self.assertEqual(fragment2.get_SMILE(), "[C]%1%2[H].[C]%1%3[H].[C]%3%4[C]%5[H].[C]%5%6[H].[C]%2%6[C]%7[O].[O]%7.[C]%4%8[O].[O]%8")



suite = unittest.TestLoader().loadTestsFromTestCase(TestFragment)
