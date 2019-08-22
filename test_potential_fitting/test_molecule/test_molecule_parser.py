import unittest
import os

from potential_fitting.molecule import molecule_parser, Atom, Fragment, Molecule
from potential_fitting.utils import SettingsReader

class TestMoleculeParser(unittest.TestCase):

    def setUpClass():
        TestMoleculeParser.monomer_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_monomer.xyz")
        TestMoleculeParser.dimer_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "NO2-_NO2-_dimer.xyz")
        TestMoleculeParser.trimer_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "NO2-_water_water_trimer.xyz")

        TestMoleculeParser.monomer_settings = SettingsReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_monomer.ini"))
        TestMoleculeParser.dimer_settings = SettingsReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_monomer.ini"))
        TestMoleculeParser.trimer_settings = SettingsReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_monomer.ini"))


    def test_xyz_to_molecules_monomer(self):
        molecules = molecule_parser.xyz_to_molecules(TestMoleculeParser.monomer_path, settings=TestMoleculeParser.monomer_settings)

        self.assertEqual(len(molecules), 2)

        for molecule in molecules:
            self.assertEqual(molecule.get_charge(), 0)
            self.assertEqual(molecule.get_spin_multiplicity(), 1)
            self.assertEqual(molecule.get_num_fragments(), 1)
            self.assertEqual(molecule.get_num_atoms(), 3)
            self.assertEqual(molecule.get_name(), "water")
            self.assertEqual(molecule.get_symmetry(), "A1B2")

        self.assertIn(Molecule([Fragment([Atom("O", "A", 1, 2, 3),
                                          Atom("H", "B", 4, 5, 6),
                                          Atom("H", "B", 7, 8, 9)
                                          ], "water", 0, 1, "O(H)H")]), molecules)

    def test_confirm_identical_results(self):
        pass

suite = unittest.TestLoader().loadTestsFromTestCase(TestMoleculeParser)
