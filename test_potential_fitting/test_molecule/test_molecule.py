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
    def test_get_fragments(self):
        molecule = Molecule([])

        # get_fragments() should return list of len 0 before any fragments added to Molecule
        self.assertEqual(len(molecule.get_fragments()), 0)

        fragment0 = Fragment([Atom("H", "A", 0, 0, 0)], "H", -3, 2, "H")

        molecule = Molecule([fragment0])

        self.assertEqual(molecule.get_fragments()[0], fragment0)
        # get_fragments() should return list of len 1 after 1 fragment added to Molecule
        self.assertEqual(len(molecule.get_fragments()), 1)
        
        fragment1 = Fragment([Atom("H", "B", 0, 0, 0), Atom("He", "C", 234, -0.1, 32), Atom("He", "D", 34, 43.53, 0)], "HHe2", 1, 1, "H[He][He]")

        molecule = Molecule([fragment0, fragment1])

        self.assertEqual(molecule.get_fragments()[0], fragment0)
        self.assertEqual(molecule.get_fragments()[1], fragment1)
        # get_fragments() should return list of len 2 after 2 fragments added to Molecule
        self.assertEqual(len(molecule.get_fragments()), 2)

    """
    Tests the get_atoms() function of the Molecule class
    """
    def test_get_atoms(self):
        molecule = Molecule([])

        # get_atoms() should return list of len 0 before any atoms added to Molecule
        self.assertEqual(len(molecule.get_atoms()), 0)

        atom0 = Atom("H", "A", 0, 0, 0)

        fragment0 = Fragment([atom0], "H", -3, 2, "H")

        molecule = Molecule([fragment0])

        self.assertEqual(molecule.get_atoms()[0], atom0)
        # get_atoms() should return list of len 1 after 1 atom added to Molecule
        self.assertEqual(len(molecule.get_atoms()), 1)

        atom1 = Atom("H", "B", 0, 0, 0)
        atom2 = Atom("He", "C", 234, -0.1, 32)
        atom3 = Atom("He", "D", 34, 43.53, 0)

        fragment1 = Fragment([atom1, atom2, atom3], "HHe2", 0, 1, "H[He][He]")

        molecule = Molecule([fragment0, fragment1])

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
        molecule = Molecule([])

        # get_charge() should return 0 before any fragments added to Molecule
        self.assertEqual(molecule.get_charge(), 0)

        molecule = Molecule([Fragment([], "Name", -1, 3, "")])

        # get_charge() should return -1 after first fragment added to Molecule
        self.assertEquals(molecule.get_charge(), -1)

        molecule = Molecule([Fragment([], "Name", -1, 3, ""), Fragment([], "Name", 3, 1, "")])

        # get_charge() should return 2 after second fragment added to Molecule
        self.assertEquals(molecule.get_charge(), 2)

    """
    Tests the get_spin_multiplicity() function of the Molecule class
    """
    def test_get_spin_multiplicity(self):
        molecule = Molecule([])

        # get_spin_multiplicity() should return 1 before any fragments added to Molecule
        self.assertEqual(molecule.get_spin_multiplicity(), 1)

        molecule = Molecule([Fragment([], "Name", -1, 3, "")])

        # get_spin_multiplicity() should return 3 after first fragment added to Molecule
        self.assertEquals(molecule.get_spin_multiplicity(), 3)

        molecule = Molecule([Fragment([], "Name", -1, 3, ""), Fragment([], "Name", 3, 1, "")])

        # get_spin_multiplicity() should return 3 after second fragment added to Molecule
        self.assertEquals(molecule.get_spin_multiplicity(), 3)

    """
    Tests the get_num_fragments() function of the Molecule class
    """
    def test_get_num_fragments(self):
        molecule = Molecule([])
        
        # get_num_fragments() should return 0 before any fragments added to molecule
        self.assertEqual(molecule.get_num_fragments(), 0)

        molecule = Molecule([Fragment([Atom("F", "A", 0, 0, 0)], "F", 1, 2, "F")])

        # get_num_fragments() should return 1 after 1 fragment added to molecule
        self.assertEqual(molecule.get_num_fragments(), 1)

        molecule = Molecule([Fragment([Atom("F", "A", 0, 0, 0)], "F", 1, 2, "F"), Fragment([Atom("F", "B", 0, 0, 3), Atom("Cl", "C", 0, 0, 0)], "FCl", 3, 2, "F[Cl]")])

        # get_num_fragments() should return 2 after 2 fragment added to molecule
        self.assertEqual(molecule.get_num_fragments(), 2)

    """
    Tests the get_num_atoms() function of the Molecule class
    """
    def test_get_num_atoms(self):
        molecule = Molecule([])
        
        # get_num_atoms() should return 0 before any atoms added to molecule
        self.assertEqual(molecule.get_num_atoms(), 0)

        molecule = Molecule([Fragment([Atom("F", "A", 0, 0, 0)], "F", 1, 2, "F")])

        # get_num_fragments() should return 1 after 1 atom added to molecule
        self.assertEqual(molecule.get_num_atoms(), 1)

        molecule = Molecule([Fragment([Atom("F", "A", 0, 0, 0)], "F", 1, 2, "F"), Fragment([Atom("F", "B", 0, 0, 3), Atom("Cl", "C", 0, 0, 0)], "FCl", 3, 2, "F[Cl]")])

        # get_num_fragments() should return 3 after 3 atoms added to molecule
        self.assertEqual(molecule.get_num_atoms(), 3)

    """
    Tests the to_xyz() function of the Molecule class
    """
    def test_to_xyz(self):

        fragment0 = Fragment([Atom("H", "A", 5, 3, 0.00343), Atom("Cl", "B", 2, 0, -13), Atom("He", "C", 6, 2, 0.343)], "HClHe", -1, 1, "H[Cl][He]")

        molecule = Molecule([fragment0])

        self.assertEqual(molecule.to_xyz(), fragment0.to_xyz()[:-1])

        fragment1 = Fragment([Atom("Ar", "D", 0.23430523424, -34, -234.5235)], "Ar", -2, 1, "[Ar]")

        molecule = Molecule([fragment0, fragment1])
        
        self.assertEqual(molecule.to_xyz(), fragment0.to_xyz() + fragment1.to_xyz()[:-1])

        fragment2 = Fragment([Atom("Xe", "E", 0, 0, 0), Atom("Br", "F", 62, 5, 0.001)], "XeBr", 0, 2, "[Xe][Br]")

        molecule = Molecule([fragment0, fragment1, fragment2])

        self.assertEqual(molecule.to_xyz(), fragment0.to_xyz() + fragment1.to_xyz() + fragment2.to_xyz()[:-1])

suite = unittest.TestLoader().loadTestsFromTestCase(TestMolecule)
