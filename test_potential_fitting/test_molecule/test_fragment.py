import unittest
import os, random

from test_potential_fitting.test_case_with_id import TestCaseWithId
from potential_fitting.molecule import Atom
from potential_fitting.molecule import Fragment
from potential_fitting.molecule import parse_training_set_file
from potential_fitting.utils import SettingsReader, Quaternion
from potential_fitting.exceptions import InvalidValueError, InconsistentValueError, XYZFormatError

"""
Test cases for Fragment class
"""
class TestFragment(TestCaseWithId):

    def test_bad_SMILE(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # tests for when a fragment is given a bad SMILE string

        # SMILE string that doesn't close all its bonds.
        with self.assertRaises(InvalidValueError):
            Fragment([Atom("O", "A", 0, 0, 0)], "frag", 0, 1, "O1")
        with self.assertRaises(InvalidValueError):
            Fragment([Atom("H", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("C", "C", 2, 2, 2)], "frag", 0, 1, "O12HC2")

        # SMILE string that indicates an atom is bonded to itself
        with self.assertRaises(InvalidValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("H", "B", 2, 2, 2)], "frag", 0, 1, "O11(H)H")

        self.test_passed = True

    def test_inconsistent_atoms_and_SMILE(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

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

        self.test_passed = True

    def test_non_positive_spin(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # test for when a fragment is given a spin multiplicity < 1

        with self.assertRaises(InvalidValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("H", "B", 2, 2, 2)], "frag", 0, 0, "O(H)H")

        with self.assertRaises(InvalidValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("H", "B", 2, 2, 2)], "frag", 0, -1, "O(H)H")

        self.test_passed = True

    def test_same_symmetry_different_type(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # test for when a fragment is given 2 atoms with the same symmetry but different atom type.

        with self.assertRaises(InconsistentValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "A", 1, 1, 1), Atom("H", "B", 2, 2, 2)], "frag", 0, 1, "O(H)H")

        with self.assertRaises(InconsistentValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("H", "A", 2, 2, 2)], "frag", 0, 1, "O(H)H")

        with self.assertRaises(InconsistentValueError):
            Fragment([Atom("O", "A", 0, 0, 0), Atom("H", "B", 1, 1, 1), Atom("C", "B", 2, 2, 2)], "frag", 0, 1, "O(H)H")

        self.test_passed = True

    def test_get_name(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "HCl", 0, 1, "")
        self.assertEqual(fragment.get_name(), "HCl")

        fragment = Fragment([], "somename!", 0, 1, "")
        self.assertEqual(fragment.get_name(), "somename!")

        self.test_passed = True

    def test_set_name(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "HCl", 0, 1, "")
        self.assertEqual(fragment.get_name(), "HCl")

        fragment.set_name("somename!")
        self.assertEqual(fragment.get_name(), "somename!")

        fragment.set_name("another name")
        self.assertEqual(fragment.get_name(), "another name")

        fragment.set_name("HCl")
        self.assertEqual(fragment.get_name(), "HCl")

        self.test_passed = True

    """
    Tests get_atoms() functions of the Fragment class
    """
    def test_get_atoms(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
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

        self.test_passed = True

    def test_get_symmetry(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "", 0, 1, "")

        # get_symmetry() should return an empty string for a fragment with no atoms
        self.assertEqual(fragment.get_symmetry(), "")

        atom0 = Atom("H", "A", 0, 0, 0)
        fragment = Fragment([atom0], "H", 0, 1, "H")

        # get_symmetry() should return "A1"
        self.assertEqual(fragment.get_symmetry(), "A1")

        atom1 = Atom("H", "A", 400, 32, 23)

        fragment = Fragment([atom0, atom1], "H2", 0, 1, "H.H")

        # get_symmetry() should return "A2"
        self.assertEqual(fragment.get_symmetry(), "A2")

        atom2 = Atom("O", "B", 34, 35, 1244)

        fragment = Fragment([atom0, atom1, atom2], "H2O", 0, 1, "H1.HO1")

        # get_symmetry() should return "A2B1"
        self.assertEqual(fragment.get_symmetry(), "A2B1")

        atom3 = Atom("O", "B", 1, 2, 3)

        fragment = Fragment([atom0, atom1, atom2, atom3], "H2O2", 0, 1, "H1.HOO1")

        # get_symmetry() should return "A2B2"
        self.assertEqual(fragment.get_symmetry(), "A2B2")

        self.test_passed = True

    def test_get_standard_symmetry(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "", 0, 1, "")

        # get_standard_symmetry() should return an empty string for a fragment with no atoms
        self.assertEqual(fragment.get_standard_symmetry(), "")

        atom0 = Atom("H", "A", 0, 0, 0)
        fragment = Fragment([atom0], "H", 0, 1, "H")

        # get_standard_symmetry() should return "A1"
        self.assertEqual(fragment.get_standard_symmetry(), "A1")

        atom1 = Atom("H", "A", 400, 32, 23)

        fragment = Fragment([atom0, atom1], "H2", 0, 1, "H.H")

        # get_standard_symmetry() should return "A2"
        self.assertEqual(fragment.get_standard_symmetry(), "A2")

        atom2 = Atom("O", "B", 34, 35, 1244)

        fragment = Fragment([atom0, atom1, atom2], "H2O", 0, 1, "H1.HO1")

        # get_standard_symmetry() should return "B1A2"
        self.assertEqual(fragment.get_standard_symmetry(), "B1A2")

        atom3 = Atom("O", "B", 1, 2, 3)

        fragment = Fragment([atom0, atom1, atom2, atom3], "H2O2", 0, 1, "H1.HOO1")

        # get_standard_symmetry() should return "B2A2"
        self.assertEqual(fragment.get_standard_symmetry(), "B2A2")

        atom4 = Atom("O", "B", 2, 3, 4)

        fragment = Fragment([atom0, atom1, atom2, atom3, atom4], "H2O3", 0, 1, "H1.HOO1O")

        # get_standard_symmetry() should return "B3A2"
        self.assertEqual(fragment.get_standard_symmetry(), "B3A2")

        self.test_passed = True

    def test_get_SMILE(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        fragment = Fragment([], "", 0, 1, "")

        # get_SMILE() should return an empty string for a fragment with no atoms
        self.assertEqual(fragment.get_SMILE(), "")

        atom0 = Atom("H", "A", 0, 0, 0)
        fragment = Fragment([atom0], "H", 0, 1, "H")

        # get_SMILE() should return "[H]"
        self.assertEqual(fragment.get_SMILE(), "[H]")

        atom1 = Atom("H", "A", 400, 32, 23)

        fragment = Fragment([atom0, atom1], "H2", 0, 1, "H.H")

        # get_SMILE() should return "[H].[H]"
        self.assertEqual(fragment.get_SMILE(), "[H].[H]")

        atom2 = Atom("O", "B", 34, 35, 1244)

        fragment = Fragment([atom0, atom1, atom2], "H2O", 0, 1, "H1.HO1")

        # get_SMILE() should return "[H]%1.[H][O]%1"
        self.assertEqual(fragment.get_SMILE(), "[H]%1.[H][O]%1")

        atom3 = Atom("O", "B", 1, 2, 3)

        fragment = Fragment([atom0, atom1, atom2, atom3], "H2O2", 0, 1, "H1.HOO1")

        # get_SMILE() should return "[H]%1.[H][O][O]%1"
        self.assertEqual(fragment.get_SMILE(), "[H]%1.[H][O][O]%1")


        # A4B2C2D4E4
        fragment1 = list(parse_training_set_file(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "bdc.xyz"), SettingsReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "bdc.ini"))))[0].get_fragments()[0]
        self.assertEqual(fragment1.get_SMILE(), "[C]%1%2[C]%3%4.[C]%5%6[C]%7%8.[C]%3%5%9.[C]%1%7%10.[H]%2.[H]%4.[H]%6.[H]%8.[C]%9%11%12.[C]%10%13%14.[O]%11.[O]%12.[O]%13.[O]%14")

        # A4B2C2D4E4
        fragment2 = list(parse_training_set_file(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "bdc2.xyz"), SettingsReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "bdc2.ini"))))[0].get_fragments()[0]
        self.assertEqual(fragment2.get_SMILE(), "[C]%1%2[H].[C]%1%3[H].[C]%3%4[C]%5[H].[C]%5%6[H].[C]%2%6[C]%7[O].[O]%7.[C]%4%8[O].[O]%8")

        self.test_passed = True

    def test_get_standard_SMILE(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        fragment = Fragment([], "", 0, 1, "")

        # get_standard_SMILE() should return an empty string for a fragment with no atoms
        self.assertEqual(fragment.get_standard_SMILE(), "")

        atom0 = Atom("H", "A", 0, 0, 0)
        fragment = Fragment([atom0], "H", 0, 1, "H")

        # get_standard_SMILE() should return "[H]"
        self.assertEqual(fragment.get_standard_SMILE(), "[H]")

        atom1 = Atom("H", "A", 400, 32, 23)

        fragment = Fragment([atom0, atom1], "H2", 0, 1, "H.H")

        # get_standard_SMILE() should return "[H].[H]"
        self.assertEqual(fragment.get_standard_SMILE(), "[H].[H]")

        atom2 = Atom("O", "B", 34, 35, 1244)

        fragment = Fragment([atom0, atom1, atom2], "H2O", 0, 1, "H1.HO1")

        # get_standard_SMILE() should return "[O]%1[H].[H]%1"
        self.assertEqual(fragment.get_standard_SMILE(), "[O]%1[H].[H]%1")

        atom3 = Atom("O", "B", 1, 2, 3)

        fragment = Fragment([atom0, atom1, atom2, atom3], "H2O2", 0, 1, "H1.HOO1")

        # get_standard_SMILE() should return "[O]%1[O]%2.[H]%1.[H]%2"
        self.assertEqual(fragment.get_standard_SMILE(), "[O]%1[O][H].[H]%1")

        self.test_passed = True


    """
    Tests the get_charge() function of the Fragment class
    """
    def test_get_charge(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "Name", 0, 2, "")
        
        self.assertEqual(fragment.get_charge(), 0)        

        fragment = Fragment([], "Name", -3, 1, "")
        
        self.assertEqual(fragment.get_charge(), -3)        
        
        fragment = Fragment([], "Name", 12, 5, "")
        
        self.assertEqual(fragment.get_charge(), 12)        

        self.test_passed = True


    """
    Tests the get_spin_multiplicity() function of the Fragment class
    """
    def test_get_spin_multiplicity(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "Name", 0, 2, "")
        
        self.assertEqual(fragment.get_spin_multiplicity(), 2)        

        fragment = Fragment([], "Name", -3, 1, "")
        
        self.assertEqual(fragment.get_spin_multiplicity(), 1)        
        
        fragment = Fragment([], "Name", 12, 5, "")
        
        self.assertEqual(fragment.get_spin_multiplicity(), 5)        

        self.test_passed = True


    """
    Tests the get_num_atoms() function of the Fragment class
    """
    def test_get_num_atoms(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
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

        self.test_passed = True

    def test_translate(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        for i in range(1000):
            start_x1 = random.random() * 100
            start_y1 = random.random() * 100
            start_z1 = random.random() * 100
            start_x2 = random.random() * 100
            start_y2 = random.random() * 100
            start_z2 = random.random() * 100

            trans_x = random.random() * 100
            trans_y = random.random() * 100
            trans_z = random.random() * 100

            fragment_ref = Fragment([Atom("H", "A", start_x1 + trans_x, start_y1 + trans_y, start_z1 + trans_z),
                                     Atom("Al", "B", start_x2 + trans_x, start_y2 + trans_y, start_z2 + trans_z)], "HAlHe", 3, 2, "H[Al]")

            fragment = Fragment([Atom("H", "A", start_x1, start_y1, start_z1),
                                      Atom("Al", "B", start_x2, start_y2, start_z2)], "HAlHe", 3, 2, "H[Al]")

            fragment.translate(trans_x, trans_y, trans_z)

            self.assertEqual(fragment, fragment_ref)

        self.test_passed = True

    def test_rotate(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([Atom("H", "A", 0, 0, 0),
                            Atom("Al", "B", 0, 0, 0)], "HAlHe", 3, 2, "H[Al]")

        for i in range(1000):

            ref_atom = Atom("H", "A", random.random() * 100, random.random() * 100, random.random() * 100)

            pre_dist0 = fragment.get_atoms()[0].distance(ref_atom)
            pre_dist1 = fragment.get_atoms()[1].distance(ref_atom)

            fragment.rotate(Quaternion.get_random_rotation_quaternion(), origin_x=ref_atom.get_x(), origin_y=ref_atom.get_y(), origin_z=ref_atom.get_z())

            post_dist0 = fragment.get_atoms()[0].distance(ref_atom)
            post_dist1 = fragment.get_atoms()[1].distance(ref_atom)

            self.assertAlmostEqual(pre_dist0, post_dist0)
            self.assertAlmostEqual(pre_dist1, post_dist1)

        self.test_passed = True

    def test_get_connectivity_matrix(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        fragment = Fragment([], "", 0, 1, "")
        self.assertEqual(fragment.get_connectivity_matrix(), [])

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "H", 0, 1, "H")
        self.assertEqual(fragment.get_connectivity_matrix(), [[False]])

        atom1 = Atom("H", "A", 400, 32, 23)

        fragment = Fragment([atom0, atom1], "H2", 0, 1, "H.H")
        self.assertEqual(fragment.get_connectivity_matrix(), [[False, False],
                                                              [False, False]])

        atom2 = Atom("O", "B", 34, 35, 1244)

        fragment = Fragment([atom0, atom1, atom2], "H2O", 0, 1, "H1.HO1")
        self.assertEqual(fragment.get_connectivity_matrix(), [[False, False, True],
                                                              [False, False, True],
                                                              [True, True, False]])

        atom3 = Atom("O", "B", 1, 2, 3)

        fragment = Fragment([atom0, atom1, atom2, atom3], "H2O2", 0, 1, "H1.HOO1")
        self.assertEqual(fragment.get_connectivity_matrix(), [[False, False, False, True],
                                                              [False, False, True, False],
                                                              [False, True, False, True],
                                                              [True, False, True, False]])

        self.test_passed = True

    def test_get_standard_connectivity_matrix(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        fragment = Fragment([], "", 0, 1, "")
        self.assertEqual(fragment.get_standard_connectivity_matrix(), [])

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "H", 0, 1, "H")
        self.assertEqual(fragment.get_standard_connectivity_matrix(), [[False]])

        atom1 = Atom("H", "A", 400, 32, 23)

        fragment = Fragment([atom0, atom1], "H2", 0, 1, "H.H")
        self.assertEqual(fragment.get_standard_connectivity_matrix(), [[False, False],
                                                                       [False, False]])

        atom2 = Atom("O", "B", 34, 35, 1244)

        fragment = Fragment([atom0, atom1, atom2], "H2O", 0, 1, "H1.HO1")
        self.assertEqual(fragment.get_standard_connectivity_matrix(), [[False, True, True],
                                                                       [True, False, False],
                                                                       [True, False, False]])

        atom3 = Atom("O", "B", 1, 2, 3)

        fragment = Fragment([atom0, atom1, atom2, atom3], "H2O2", 0, 1, "H1.HOO1")
        self.assertEqual(fragment.get_standard_connectivity_matrix(), [[False, True, False, True],
                                                                       [True, False, True, False],
                                                                       [False, True, False, False],
                                                                       [True, False, False, False]])

        self.test_passed = True

    def test_get_excluded_pairs(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        fragment = Fragment([], "", 0, 1, "")
        self.assertEqual(fragment.get_excluded_pairs(), [[],
                                                         [],
                                                         []])

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "H", 0, 1, "H")
        self.assertEqual(fragment.get_excluded_pairs(), [[],
                                                         [],
                                                         []])

        atom1 = Atom("H", "A", 400, 32, 23)

        fragment = Fragment([atom0, atom1], "H2", 0, 1, "H.H")
        self.assertEqual(fragment.get_excluded_pairs(), [[],
                                                         [],
                                                         []])

        atom2 = Atom("O", "B", 34, 35, 1244)

        fragment = Fragment([atom0, atom1, atom2], "H2O", 0, 1, "H1.HO1")
        excluded12, excluded13, excluded14 = fragment.get_excluded_pairs()
        self.assertEqual(len(excluded12), 2)
        self.assertEqual(len(excluded13), 1)
        self.assertEqual(len(excluded14), 0)
        self.assertIn([0, 2], excluded12)
        self.assertIn([1, 2], excluded12)
        self.assertIn([0, 1], excluded13)

        atom3 = Atom("O", "B", 1, 2, 3)

        fragment = Fragment([atom0, atom1, atom2, atom3], "H2O2", 0, 1, "H1.HOO1")
        excluded12, excluded13, excluded14 = fragment.get_excluded_pairs()
        self.assertEqual(len(excluded12), 3)
        self.assertEqual(len(excluded13), 2)
        self.assertEqual(len(excluded14), 1)
        self.assertIn([0, 3], excluded12)
        self.assertIn([1, 2], excluded12)
        self.assertIn([2, 3], excluded12)
        self.assertIn([0, 2], excluded13)
        self.assertIn([1, 3], excluded13)
        self.assertIn([0, 1], excluded14)

        self.test_passed = True

    """
    Test the to_xyz() function of the Fragment class
    """
    def test_to_xyz(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "fragment", -2, 2, "")

        # to_xyz() should return empty string when no atoms are added to fragment
        self.assertEqual(fragment.to_xyz(), "")

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "fragment", -2, 2, "H")

        # to_xyz() should return string of first atom after only 1 atom added
        self.assertEqual(fragment.to_xyz(), atom0.to_xyz() + "\n")

        atom1 = Atom("Cl", "B", 5, 7, -3)

        fragment = Fragment([atom0, atom1], "fragment", -2, 2, "H[Cl]")

        # to_xyz() should return string of 2 atoms after 2nd atom added
        self.assertEqual(fragment.to_xyz(), atom0.to_xyz() + "\n" + atom1.to_xyz() + "\n")

        atom2 = Atom("Xe", "C", 10.234235, -0.00000234, 2.353523)

        fragment = Fragment([atom0, atom1, atom2], "fragment", -2, 2, "H[Cl][Xe]")

        # to_xyz() should return string of 3 atoms after only 3rd atom added
        self.assertEqual(fragment.to_xyz(), atom0.to_xyz() + "\n" + atom1.to_xyz() + "\n" + atom2.to_xyz() + "\n")

        self.test_passed = True

    def test_to_standard_xyz(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "fragment", -2, 2, "")

        # to_xyz() should return empty string when no atoms are added to fragment
        self.assertEqual(fragment.to_standard_xyz(), "")

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "fragment", -2, 2, "H")

        # to_standard_xyz() should return string of first atom after only 1 atom added
        self.assertEqual(fragment.to_standard_xyz(), atom0.to_xyz() + "\n")

        atom1 = Atom("Cl", "B", 5, 7, -3)

        fragment = Fragment([atom0, atom1], "fragment", -2, 2, "H[Cl]")

        # to_standard_xyz() should return string of 2 atoms after 2nd atom added
        self.assertEqual(fragment.to_standard_xyz(), atom1.to_xyz() + "\n" + atom0.to_xyz() + "\n")

        atom2 = Atom("Xe", "C", 10.234235, -0.00000234, 2.353523)

        fragment = Fragment([atom0, atom1, atom2], "fragment", -2, 2, "H[Cl][Xe]")

        # to_standard_xyz() should return string of 3 atoms after only 3rd atom added
        self.assertEqual(fragment.to_standard_xyz(), atom2.to_xyz() + "\n" + atom1.to_xyz() + "\n" + atom0.to_xyz() + "\n")

        self.test_passed = True

    """
    Test the to_xyz() function of the Fragment class
    """
    def test_to_ghost_xyz(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "fragment", -2, 2, "")

        # to_ghost_xyz() should return empty string when no atoms are added to fragment
        self.assertEqual(fragment.to_ghost_xyz(), "")

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "fragment", -2, 2, "H")

        # to_ghost_xyz() should return string of first atom after only 1 atom added
        self.assertEqual(fragment.to_ghost_xyz(), atom0.to_ghost_xyz() + "\n")

        atom1 = Atom("Cl", "B", 5, 7, -3)

        fragment = Fragment([atom0, atom1], "fragment", -2, 2, "H[Cl]")

        # to_ghost_xyz() should return string of 2 atoms after 2nd atom added
        self.assertEqual(fragment.to_ghost_xyz(), atom0.to_ghost_xyz() + "\n" + atom1.to_ghost_xyz() + "\n")

        atom2 = Atom("Xe", "C", 10.234235, -0.00000234, 2.353523)

        fragment = Fragment([atom0, atom1, atom2], "fragment", -2, 2, "H[Cl][Xe]")

        # to_ghost_xyz() should return string of 3 atoms after only 3rd atom added
        self.assertEqual(fragment.to_ghost_xyz(), atom0.to_ghost_xyz() + "\n" + atom1.to_ghost_xyz() + "\n" + atom2.to_ghost_xyz() + "\n")

        self.test_passed = True

    def test_to_ghost_standard_xyz(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "fragment", -2, 2, "")

        # to_standard_ghost_xyz() should return empty string when no atoms are added to fragment
        self.assertEqual(fragment.to_standard_ghost_xyz(), "")

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "fragment", -2, 2, "H")

        # to_standard_ghost_xyz() should return string of first atom after only 1 atom added
        self.assertEqual(fragment.to_standard_ghost_xyz(), atom0.to_ghost_xyz() + "\n")

        atom1 = Atom("Cl", "B", 5, 7, -3)

        fragment = Fragment([atom0, atom1], "fragment", -2, 2, "H[Cl]")

        # to_standard_ghost_xyz() should return string of 2 atoms after 2nd atom added
        self.assertEqual(fragment.to_standard_ghost_xyz(), atom1.to_ghost_xyz() + "\n" + atom0.to_ghost_xyz() + "\n")

        atom2 = Atom("Xe", "C", 10.234235, -0.00000234, 2.353523)

        fragment = Fragment([atom0, atom1, atom2], "fragment", -2, 2, "H[Cl][Xe]")

        # to_standard_ghost_xyz() should return string of 3 atoms after only 3rd atom added
        self.assertEqual(fragment.to_standard_ghost_xyz(), atom2.to_ghost_xyz() + "\n" + atom1.to_ghost_xyz() + "\n" + atom0.to_ghost_xyz() + "\n")

        self.test_passed = True

    def test_read_xyz(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment.read_xyz("H 0 0 0", "fragment", -2, 2, "H", "A1")
        fragment_ref = Fragment([atom0], "fragment", -2, 2, "H")

        self.assertEqual(fragment, fragment_ref)

        atom1 = Atom("Cl", "B", 5, 7, -3)

        fragment = Fragment.read_xyz("H 0 0 0\nCl 5 7 -3", "fragment", -2, 2, "H[Cl]", "A1B1")
        fragment_ref = Fragment([atom0, atom1], "fragment", -2, 2, "H[Cl]")

        self.assertEqual(fragment, fragment_ref)

        atom2 = Atom("Xe", "C", 10.234235, -0.00000234, 2.353523)

        fragment = Fragment.read_xyz("H 0 0 0\nCl 5 7 -3\nXe 10.234235 -0.00000234 2.353523", "fragment", -2, 2, "H[Cl][Xe]", "A1B1C1")
        fragment_ref = Fragment([atom0, atom1, atom2], "fragment", -2, 2, "H[Cl][Xe]")

        self.assertEqual(fragment, fragment_ref)

        self.test_passed = True

    def test_read_xyz_bad_line_format(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # too few arguments on atom line

        with self.assertRaises(XYZFormatError):
            Fragment.read_xyz("H 0 0", "fragment", -2, 2, "H", "A1")

        # too many arguments on atom line

        with self.assertRaises(XYZFormatError):
            Fragment.read_xyz("H 0 0 0 0", "fragment", -2, 2, "H", "A1")

        self.test_passed = True

    def test_read_xyz_bad_symmetry(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        # more atoms in symmetry than string

        with self.assertRaises(InconsistentValueError):
            Fragment.read_xyz("H 0 0 0", "fragment", -2, 2, "H", "A2")
        with self.assertRaises(InconsistentValueError):
            Fragment.read_xyz("H 0 0 0", "fragment", -2, 2, "H", "A1B1")

        # fewer atoms in symmetry than string

        with self.assertRaises(InconsistentValueError):
            Fragment.read_xyz("H 0 0 0\nH 0 0 0", "fragment", -2, 2, "H", "A1")
        with self.assertRaises(InconsistentValueError):
            Fragment.read_xyz("H 0 0 0\n O 0 0 0\n O 0 0 0", "fragment", -2, 2, "H", "A1B1")

        self.test_passed = True

    def test_get_standard_order(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "fragment", -2, 2, "")
        self.assertEqual(fragment.get_standard_order(), [])

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "fragment", -2, 2, "H")
        self.assertEqual(fragment.get_standard_order(), [atom0])

        atom1 = Atom("Cl", "B", 5, 7, -3)

        fragment = Fragment([atom0, atom1], "fragment", -2, 2, "H[Cl]")
        self.assertEqual(fragment.get_standard_order(), [atom1, atom0])

        atom2 = Atom("Xe", "C", 10.234235, -0.00000234, 2.353523)

        fragment = Fragment([atom0, atom1, atom2], "fragment", -2, 2, "H[Cl][Xe]")
        self.assertEqual(fragment.get_standard_order(), [atom2, atom1, atom0])

        self.test_passed = True

    def test_confirm_standard_order(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "fragment", -2, 2, "")
        self.assertTrue(fragment.confirm_standard_order())

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "fragment", -2, 2, "H")
        self.assertTrue(fragment.confirm_standard_order())

        atom1 = Atom("Cl", "B", 5, 7, -3)

        fragment = Fragment([atom0, atom1], "fragment", -2, 2, "H[Cl]")
        self.assertFalse(fragment.confirm_standard_order())

        fragment = Fragment([atom1, atom0], "fragment", -2, 2, "[Cl]H")
        self.assertTrue(fragment.confirm_standard_order())

        atom2 = Atom("Xe", "C", 10.234235, -0.00000234, 2.353523)

        fragment = Fragment([atom0, atom1, atom2], "fragment", -2, 2, "H[Cl][Xe]")
        self.assertFalse(fragment.confirm_standard_order())
        fragment = Fragment([atom0, atom2, atom1], "fragment", -2, 2, "H[Xe][Cl]")
        self.assertFalse(fragment.confirm_standard_order())
        fragment = Fragment([atom1, atom0, atom2], "fragment", -2, 2, "[Cl]H[Xe]")
        self.assertFalse(fragment.confirm_standard_order())
        fragment = Fragment([atom1, atom2, atom0], "fragment", -2, 2, "[Cl][Xe]H")
        self.assertFalse(fragment.confirm_standard_order())
        fragment = Fragment([atom2, atom0, atom1], "fragment", -2, 2, "[Xe]H[Cl]")
        self.assertFalse(fragment.confirm_standard_order())

        fragment = Fragment([atom2, atom1, atom0], "fragment", -2, 2, "[Xe][Cl]H")
        self.assertTrue(fragment.confirm_standard_order())

        self.test_passed = True

    def test_confirm_symmetry_class(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        fragment = Fragment([], "fragment", -2, 2, "")
        self.assertTrue(fragment.confirm_symmetry_class()[0])
        self.assertEqual(fragment.confirm_symmetry_class()[1], "")
        self.assertEqual(fragment.confirm_symmetry_class()[2], "")

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment = Fragment([atom0], "fragment", -2, 2, "H")
        self.assertTrue(fragment.confirm_symmetry_class()[0])
        self.assertEqual(fragment.confirm_symmetry_class()[1], "A1")
        self.assertEqual(fragment.confirm_symmetry_class()[2], "A1")

        atom1 = Atom("O", "B", 1, 1, 1)

        fragment = Fragment([atom1], "fragment", -2, 2, "O")
        self.assertTrue(fragment.confirm_symmetry_class()[0])
        self.assertEqual(fragment.confirm_symmetry_class()[1], "A1")
        self.assertEqual(fragment.confirm_symmetry_class()[2], "A1")

        atom2 = Atom("H", "A", 2, 2, 2)

        fragment = Fragment([atom0, atom2, atom1], "fragment", -2, 2, "H1.HO1")
        self.assertTrue(fragment.confirm_symmetry_class()[0])
        self.assertEqual(fragment.confirm_symmetry_class()[1], "A1B2")
        self.assertEqual(fragment.confirm_symmetry_class()[2], "A1B2")

        fragment = Fragment([atom0, atom2, atom1], "fragment", -2, 2, "HHO")
        self.assertFalse(fragment.confirm_symmetry_class()[0])
        self.assertEqual(fragment.confirm_symmetry_class()[1], "A1B1C1")
        self.assertEqual(fragment.confirm_symmetry_class()[2], "A1B2")

        atom2 = Atom("H", "C", 2, 2, 2)

        fragment = Fragment([atom0, atom2, atom1], "fragment", -2, 2, "H1.HO1")
        self.assertFalse(fragment.confirm_symmetry_class()[0])
        self.assertEqual(fragment.confirm_symmetry_class()[1], "A1B2")
        self.assertEqual(fragment.confirm_symmetry_class()[2], "A1B1C1")

        self.test_passed = True

    def test_get_standard_copy(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        atom_std1 = Atom("O", "A", 1, 1, 1)
        atom_std2 = Atom("C", "B", 2, 2, 2)
        atom_std3 = Atom("H", "C", 0, 0, 0)

        standard_frag = Fragment([atom_std1, atom_std2, atom_std3], "fragment", 0, 1, "O(C)H")

        self.assertEqual(standard_frag.get_standard_copy(), standard_frag)

        atom0 = Atom("H", "A", 0, 0, 0)
        atom1 = Atom("C", "B", 2, 2, 2)
        atom2 = Atom("O", "C", 1, 1, 1)

        frag = Fragment([atom0, atom1, atom2], "fragment", 0, 1, "H1.CO1")

        self.assertEqual(frag.get_standard_copy(), standard_frag)

        atom0 = Atom("C", "A", 2, 2, 2)
        atom1 = Atom("H", "B", 0, 0, 0)
        atom2 = Atom("O", "C", 1, 1, 1)

        frag = Fragment([atom0, atom1, atom2], "fragment", 0, 1, "C1.HO1")

        self.assertEqual(frag.get_standard_copy(), standard_frag)

        self.test_passed = True

    def test_get_reorder_copy(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        atom_std1 = Atom("O", "A", 1, 1, 1)
        atom_std2 = Atom("C", "B", 2, 2, 2)
        atom_std3 = Atom("H", "C", 0, 0, 0)

        standard_frag = Fragment([atom_std1, atom_std2, atom_std3], "fragment", 0, 1, "O(C)H")

        self.assertEqual(standard_frag.get_reorder_copy("O(C)H"), standard_frag)

        atom0 = Atom("H", "A", 0, 0, 0)
        atom1 = Atom("C", "B", 2, 2, 2)
        atom2 = Atom("O", "C", 1, 1, 1)

        frag = Fragment([atom0, atom1, atom2], "fragment", 0, 1, "H1.CO1")

        self.assertEqual(standard_frag.get_reorder_copy("H1.CO1"), frag)

        atom0 = Atom("H", "A", 0, 0, 0)
        atom1 = Atom("O", "B", 1, 1, 1)
        atom2 = Atom("C", "C", 2, 2, 2)

        frag = Fragment([atom0, atom1, atom2], "fragment", 0, 1, "HOC")

        self.assertEqual(standard_frag.get_reorder_copy("HOC"), frag)

        self.test_passed = True

    def test_eq_and_ne(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        atom0 = Atom("O", "A", 1, 1, 1)

        atom1 = Atom("H", "B", 0, 0, 0)

        atom2 = Atom("H", "B", 2, 2, 2)

        frag0 = Fragment([atom0, atom1, atom2], "fragment", 0, 1, "O(H)H")
        frag1 = Fragment([atom0, atom1, atom2], "fragment", 0, 1, "O(H)H")

        self.assertEqual(frag0, frag1)

        frag2 = Fragment([atom0, atom1, atom2], "fragment", 0, 2, "O(H)H")

        self.assertNotEqual(frag0, frag2)

        frag3 = Fragment([atom0, atom1, atom2], "fragment", 1, 1, "O(H)H")

        self.assertNotEqual(frag0, frag3)

        atom3 = Atom("H", "C", 2, 2, 2)

        frag4 = Fragment([atom0, atom1, atom3], "fragment", 0, 1, "O(H)H")

        self.assertNotEqual(frag0, frag4)

        self.test_passed = True



suite = unittest.TestLoader().loadTestsFromTestCase(TestFragment)
