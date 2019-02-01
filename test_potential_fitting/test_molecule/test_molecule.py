import unittest

from potential_fitting.molecule import Atom
from potential_fitting.molecule import Fragment
from potential_fitting.molecule import Molecule

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
        atom0 = Atom("H", "A", 0, 0, 0)

        fragment0.add_atom(atom0)

        molecule.add_fragment(fragment0)

        self.assertEqual(molecule.get_fragments()[0], fragment0)
        # get_fragments() should return list of len 1 after 1 fragment added to Molecule
        self.assertEqual(len(molecule.get_fragments()), 1)
        
        fragment1 = Fragment("HHe2", 1, 1)
        atom1 = Atom("H", "B", 0, 0, 0)
        atom2 = Atom("He", "C", 234, -0.1, 32)
        atom3 = Atom("He", "D", 34, 43.53, 0)

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
        atom0 = Atom("H", "A", 0, 0, 0)

        fragment0.add_atom(atom0)

        molecule.add_fragment(fragment0)

        self.assertEqual(molecule.get_atoms()[0], atom0)
        # get_atoms() should return list of len 1 after 1 atom added to Molecule
        self.assertEqual(len(molecule.get_atoms()), 1)
        
        
        fragment1 = Fragment("HHe2", 0, 1)
        atom1 = Atom("H", "B", 0, 0, 0)
        atom2 = Atom("He", "C", 234, -0.1, 32)
        atom3 = Atom("He", "D", 34, 43.53, 0)

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
        atom0 = Atom("F", "A", 0, 0, 0)

        fragment0.add_atom(atom0)

        molecule.add_fragment(fragment0)

        # get_num_fragments() should return 1 after 1 fragment added to molecule
        self.assertEqual(molecule.get_num_fragments(), 1)

        fragment1 = Fragment("FCl", 3, 2)
        atom1 = Atom("F", "B", 0, 0, 3)
        atom2 = Atom("Cl", "C", 0, 0, 0)

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
        atom0 = Atom("F", "A", 0, 0, 0)

        fragment0.add_atom(atom0)

        molecule.add_fragment(fragment0)

        # get_num_fragments() should return 1 after 1 atom added to molecule
        self.assertEqual(molecule.get_num_atoms(), 1)

        fragment1 = Fragment("FCl", -1, 2)
        atom1 = Atom("F", "B", 0, 0, 3)
        atom2 = Atom("Cl", "C", 0, 0, 0)

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
        fragment0.add_atom(Atom("H", "A", 5, 3, 0.00343))
        fragment0.add_atom(Atom("Cl", "B", 2, 0, -13))
        fragment0.add_atom(Atom("He", "C", 6, 2, 0.343))

        molecule.add_fragment(fragment0)

        self.assertEqual(molecule.to_xyz(), fragment0.to_xyz()[:-1])

        fragment1 = Fragment("Ar", -2, 1)
        fragment1.add_atom(Atom("Ar", "D", 0.23430523424, -34, -234.5235))

        molecule.add_fragment(fragment1)
        
        self.assertEqual(molecule.to_xyz(), fragment0.to_xyz() + fragment1.to_xyz()[:-1])

        fragment2 = Fragment("XeBr", 0, 2)
        fragment2.add_atom(Atom("Xe", "E", 0, 0, 0))
        fragment2.add_atom(Atom("Br", "F", 62, 5, 0.001))

        molecule.add_fragment(fragment2)
        self.assertEqual(molecule.to_xyz(), fragment0.to_xyz() + fragment1.to_xyz() + fragment2.to_xyz()[:-1])

suite = unittest.TestLoader().loadTestsFromTestCase(TestMolecule)
