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
        atom = Atom("H", 0, 0, 0)
        self.assertEqual(atom.get_name(), "H")
        
        atom = Atom("Cl", 0, 0, 0)
        self.assertEqual(atom.get_name(), "Cl")
        
        atom = Atom("He", 0, 0, 0)
        self.assertEqual(atom.get_name(), "He")
        
        atom = Atom("Ar", 0, 0, 0)
        self.assertEqual(atom.get_name(), "Ar")

    """
    Tests the to_xyz() function of the Atom class
    """ 
    def test_to_xyz(self):
        atom = Atom("H", 0, 0, 0)
        self.assertEqual(atom.to_xyz(), "H    0.00000000000000e+00   0.00000000000000e+00   0.00000000000000e+00");
        
        atom = Atom("Ar", 0.34234, -0.342, 1.2334)
        self.assertEqual(atom.to_xyz(), "Ar   3.42340000000000e-01  -3.42000000000000e-01   1.23340000000000e+00");
        
        atom = Atom("He", 12.34, 105.34, -0.00432)
        self.assertEqual(atom.to_xyz(), "He   1.23400000000000e+01   1.05340000000000e+02  -4.32000000000000e-03");
        
        atom = Atom("C", -2523452, 0.0003453, 34534)
        self.assertEqual(atom.to_xyz(), "C   -2.52345200000000e+06   3.45300000000000e-04   3.45340000000000e+04");

        atom = Atom("Cl", 0.00000000000068, 0.74576456847823583, 34534262462472457756745)
        self.assertEqual(atom.to_xyz(), "Cl   6.80000000000000e-13   7.45764568478236e-01   3.45342624624725e+22");

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

        atom0 = Atom("H", 0, 0, 0)

        fragment.add_atom(atom0)

        self.assertEqual(fragment.get_atoms()[0], atom0)
        # get_atom() should return list of length 1 after 1 atom added to fragment
        self.assertEqual(len(fragment.get_atoms()), 1)

        atom1 = Atom("Cl", 400, 32, 23)    
    
        fragment.add_atom(atom1)
        
        self.assertEqual(fragment.get_atoms()[0], atom1)
        self.assertEqual(fragment.get_atoms()[1], atom0)
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
        
        fragment.add_atom(Atom("H", 0, 0, 0))

        # get_num_atoms() should return 1 after single atom added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 1)
        
        fragment.add_atom(Atom("Al", 0, 0, 0))

        # get_num_atoms() should return 2 after second atom added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 2)
        
        fragment.add_atom(Atom("He", 0, 0, 0))

        # get_num_atoms() should return 3 after third atom added to Fragment
        self.assertEqual(fragment.get_num_atoms(), 3)

    """
    Test the to_xyz() function of the Fragment class
    """
    def test_to_xyz(self):
        fragment = Fragment("HClXe", -2, 2);

        # to_xyz() should return empty string when no atoms are added to fragment
        self.assertEqual(fragment.to_xyz(), "");
        
        atom0 = Atom("H", 0, 0, 0)

        fragment.add_atom(atom0)

        # to_xyz() should return string of first atom after only 1 atom added
        self.assertEqual(fragment.to_xyz(), atom0.to_xyz() + "\n")
        
        atom1 = Atom("Cl", 5, 7, -3)

        fragment.add_atom(atom1)

        # to_xyz() should return string of 2 atoms after 2nd atom added
        self.assertEqual(fragment.to_xyz(), atom1.to_xyz() + "\n" + atom0.to_xyz() + "\n")
    
        atom2 = Atom("Xe", 10.234235, -0.00000234, 2.353523)

        fragment.add_atom(atom2)

        # to_xyz() should return string of 3 atoms after only 3rd atom added
        self.assertEqual(fragment.to_xyz(), atom1.to_xyz() + "\n" + atom0.to_xyz() + "\n" + atom2.to_xyz() + "\n")

"""
Test cases for molecule class
"""
class TestMolecule(unittest.TestCase):
    """
    Test the add_fragment() and get_fragments() functions of the Molecule class
    """
    def test_add_fragment_and_get_fragments(self):
        molecule = Molecule()

        # get_fragments() should return list of len 0 before any fragments added to Molecule
        self.assertEqual(len(molecule.get_fragments()), 0)

        fragment0 = Fragment("H", -3, 2)
        atom0 = Atom("H", 0, 0, 0)

        fragment0.add_atom(atom0)

        molecule.add_fragment(fragment0)

        self.assertEqual(molecule.get_fragments()[0], fragment0)
        # get_fragments() should return list of len 1 after 1 fragment added to Molecule
        self.assertEqual(len(molecule.get_fragments()), 1)
        
        fragment1 = Fragment("HHe2", 1, 1)
        atom1 = Atom("H", 0, 0, 0)
        atom2 = Atom("He", 234, -0.1, 32)
        atom3 = Atom("He", 34, 43.53, 0)

        fragment1.add_atom(atom1)
        fragment1.add_atom(atom2)
        fragment1.add_atom(atom3)

        molecule.add_fragment(fragment1)

        self.assertEqual(molecule.get_fragments()[0], fragment0)
        self.assertEqual(molecule.get_fragments()[1], fragment1)
        # get_fragments() should return list of len 2 after 2 fragments added to Molecule
        self.assertEqual(len(molecule.get_fragments()), 2)

    """
    Tests the get_atoms() function of the Molecule class
    """
    def test_get_atoms(self):
        molecule = Molecule()

        # get_atoms() should return list of len 0 before any atoms added to Molecule
        self.assertEqual(len(molecule.get_atoms()), 0)

        fragment0 = Fragment("H", 0, 1)
        atom0 = Atom("H", 0, 0, 0)

        fragment0.add_atom(atom0)

        molecule.add_fragment(fragment0)

        self.assertEqual(molecule.get_atoms()[0], atom0)
        # get_atoms() should return list of len 1 after 1 atom added to Molecule
        self.assertEqual(len(molecule.get_atoms()), 1)
        
        
        fragment1 = Fragment("HHe2", 0, 1)
        atom1 = Atom("H", 0, 0, 0)
        atom2 = Atom("He", 234, -0.1, 32)
        atom3 = Atom("He", 34, 43.53, 0)

        fragment1.add_atom(atom1)
        fragment1.add_atom(atom2)
        fragment1.add_atom(atom3)

        molecule.add_fragment(fragment1)

        self.assertEqual(molecule.get_atoms()[0], atom0)
        self.assertEqual(molecule.get_atoms()[1], atom1)
        self.assertEqual(molecule.get_atoms()[2], atom2)
        self.assertEqual(molecule.get_atoms()[3], atom3)
        # get_atoms() should return list of len 4 after 4 atoms added to molecule
        self.assertEqual(len(molecule.get_atoms()), 4)

    """
    Tests the get_charge() function of the Molecule class
    """
    def test_get_charge(self):
        molecule = Molecule()

        # get_charge() should return 0 before any fragments added to Molecule
        self.assertEqual(molecule.get_charge(), 0)
        
        fragment0 = Fragment("Name", -1, 3)

        molecule.add_fragment(fragment0)

        # get_charge() should return -1 after first fragment added to Molecule
        self.assertEquals(molecule.get_charge(), -1)

        fragment1 = Fragment("Name", 3, 1)

        molecule.add_fragment(fragment1)

        # get_charge() should return 2 after second fragment added to Molecule
        self.assertEquals(molecule.get_charge(), 2)

    """
    Tests the get_spin_multiplicity() function of the Molecule class
    """
    def test_get_spin_multiplicity(self):
        molecule = Molecule()

        # get_spin_multiplicity() should return 1 before any fragments added to Molecule
        self.assertEqual(molecule.get_spin_multiplicity(), 1)
        
        fragment0 = Fragment("Name", -1, 3)

        molecule.add_fragment(fragment0)

        # get_spin_multiplicity() should return 3 after first fragment added to Molecule
        self.assertEquals(molecule.get_spin_multiplicity(), 3)

        fragment1 = Fragment("Name", 3, 1)

        molecule.add_fragment(fragment1)

        # get_spin_multiplicity() should return 3 after second fragment added to Molecule
        self.assertEquals(molecule.get_spin_multiplicity(), 3)

    """
    Tests the get_num_fragments() function of the Molecule class
    """
    def test_get_num_fragments(self):
        molecule = Molecule()
        
        # get_num_fragments() should return 0 before any fragments added to molecule
        self.assertEqual(molecule.get_num_fragments(), 0)
        
        fragment0 = Fragment("F", 1, 2)
        atom0 = Atom("F", 0, 0, 0)

        fragment0.add_atom(atom0)

        molecule.add_fragment(fragment0)

        # get_num_fragments() should return 1 after 1 fragment added to molecule
        self.assertEqual(molecule.get_num_fragments(), 1)

        fragment1 = Fragment("FCl", 3, 2)
        atom1 = Atom("F", 0, 0, 3)
        atom2 = Atom("Cl", 0, 0, 0)

        fragment1.add_atom(atom1)
        fragment1.add_atom(atom2)

        molecule.add_fragment(fragment1)

        # get_num_fragments() should return 2 after 2 fragment added to molecule
        self.assertEqual(molecule.get_num_fragments(), 2)

    """
    Tests the get_num_atoms() function of the Molecule class
    """
    def test_get_num_atoms(self):
        molecule = Molecule()
        
        # get_num_atoms() should return 0 before any atoms added to molecule
        self.assertEqual(molecule.get_num_atoms(), 0)
        
        fragment0 = Fragment("F", -1, 2)
        atom0 = Atom("F", 0, 0, 0)

        fragment0.add_atom(atom0)

        molecule.add_fragment(fragment0)

        # get_num_fragments() should return 1 after 1 atom added to molecule
        self.assertEqual(molecule.get_num_atoms(), 1)

        fragment1 = Fragment("FCl", -1, 2)
        atom1 = Atom("F", 0, 0, 3)
        atom2 = Atom("Cl", 0, 0, 0)

        fragment1.add_atom(atom1)
        fragment1.add_atom(atom2)

        molecule.add_fragment(fragment1)

        # get_num_fragments() should return 3 after 3 atoms added to molecule
        self.assertEqual(molecule.get_num_atoms(), 3)

    """
    Tests the to_xyz() function of the Molecule class
    """
    def test_to_xyz(self):
        molecule = Molecule()

        fragment0 = Fragment("HClHe", -1, 1)
        fragment0.add_atom(Atom("H", 5, 3, 0.00343))
        fragment0.add_atom(Atom("Cl", 2, 0, -13))
        fragment0.add_atom(Atom("He", 6, 2, 0.343))

        molecule.add_fragment(fragment0)

        self.assertEqual(molecule.to_xyz(), fragment0.to_xyz()[:-1])

        fragment1 = Fragment("Ar", -2, 1)
        fragment1.add_atom(Atom("Ar", 0.23430523424, -34, -234.5235))

        molecule.add_fragment(fragment1)
        
        self.assertEqual(molecule.to_xyz(), fragment1.to_xyz() + fragment0.to_xyz()[:-1])

        fragment2 = Fragment("XeBr", 0, 2)
        fragment2.add_atom(Atom("Xe", 0, 0, 0))
        fragment2.add_atom(Atom("Br", 62, 5, 0.001))

        molecule.add_fragment(fragment2)
        self.assertEqual(molecule.to_xyz(), fragment1.to_xyz() + fragment0.to_xyz() + fragment2.to_xyz()[:-1])

atomSuite = unittest.TestLoader().loadTestsFromTestCase(TestAtom)
fragmentSuite = unittest.TestLoader().loadTestsFromTestCase(TestFragment)
moleculeSuite = unittest.TestLoader().loadTestsFromTestCase(TestMolecule)

suite = unittest.TestSuite([atomSuite, fragmentSuite, moleculeSuite])
