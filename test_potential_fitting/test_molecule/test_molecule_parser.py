import unittest
import os

from test_potential_fitting.test_case_with_id import TestCaseWithId
from potential_fitting.molecule import molecule_parser, Atom, Fragment, Molecule
from potential_fitting.utils import SettingsReader

class TestMoleculeParser(TestCaseWithId):

    def setUpClass():
        TestMoleculeParser.monomer_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_monomer.xyz")
        TestMoleculeParser.dimer_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "NO2-_NO2-_dimer.xyz")
        TestMoleculeParser.trimer_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "NO2-_water_water_trimer.xyz")

        TestMoleculeParser.monomer_settings = SettingsReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "water_monomer.ini"))
        TestMoleculeParser.dimer_settings = SettingsReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "NO2-_NO2-_dimer.ini"))
        TestMoleculeParser.trimer_settings = SettingsReader(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "NO2-_water_water_trimer.ini"))


    def test_xyz_to_molecules_monomer(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
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

        self.assertIn(Molecule([Fragment([Atom("O", "A", 2, 3, 4),
                                          Atom("H", "B", 5, 6, 7),
                                          Atom("H", "B", 8, 9, 10)
                                          ], "water", 0, 1, "O(H)H")]), molecules)

        self.test_passed = True

    def test_xyz_to_molecules_dimer(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        molecules = molecule_parser.xyz_to_molecules(TestMoleculeParser.dimer_path, settings=TestMoleculeParser.dimer_settings)

        self.assertEqual(len(molecules), 2)

        for molecule in molecules:
            self.assertEqual(molecule.get_charge(), -2)
            self.assertEqual(molecule.get_spin_multiplicity(), 1)
            self.assertEqual(molecule.get_num_fragments(), 2)
            self.assertEqual(molecule.get_num_atoms(), 6)
            self.assertEqual(molecule.get_name(), "NO2--NO2-")
            self.assertEqual(molecule.get_symmetry(), "A1B2_A1B2")

        self.assertIn(Molecule([Fragment([Atom("N", "A", 1, 2, 3),
                                          Atom("O", "B", 4, 5, 6),
                                          Atom("O", "B", 7, 8, 9),
                                          ], "NO2-", -1, 1, "N(O)O"),
                                Fragment([Atom("N", "A", 2, 3, 4),
                                          Atom("O", "B", 5, 6, 7),
                                          Atom("O", "B", 8, 9, 10)
                                          ], "NO2-", -1, 1, "N(O)O")
                                ]), molecules)

        self.assertIn(Molecule([Fragment([Atom("N", "A", -1, -2, -3.6),
                                          Atom("O", "B", -4, -5, -6e2),
                                          Atom("O", "B", -7, -8, -9),
                                          ], "NO2-", -1, 1, "N(O)O"),
                                Fragment([Atom("N", "A", -2, -3, -4),
                                          Atom("O", "B", -5, -6, -7),
                                          Atom("O", "B", -8, -9, -10)
                                          ], "NO2-", -1, 1, "N(O)O")
                                ]), molecules)

        self.test_passed = True

    def test_xyz_to_molecules_trimer(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        molecules = molecule_parser.xyz_to_molecules(TestMoleculeParser.trimer_path, settings=TestMoleculeParser.trimer_settings)

        self.assertEqual(len(molecules), 1)

        for molecule in molecules:
            self.assertEqual(molecule.get_charge(), -1)
            self.assertEqual(molecule.get_spin_multiplicity(), 1)
            self.assertEqual(molecule.get_num_fragments(), 3)
            self.assertEqual(molecule.get_num_atoms(), 9)
            self.assertEqual(molecule.get_name(), "NO2--water-water")
            self.assertEqual(molecule.get_symmetry(), "A1B2_C1D2_C1D2")

        self.assertIn(Molecule([Fragment([Atom("N", "A", 1, 2, 3),
                                          Atom("O", "B", 4, 5, 6),
                                          Atom("O", "B", 7, 8, 9),
                                          ], "NO2-", -1, 1, "N(O)O"),
                                Fragment([Atom("O", "C", 2, 3, 4),
                                          Atom("H", "D", 5, 6, 7),
                                          Atom("H", "D", 8, 9, 1)
                                          ], "water", 0, 1, "O(H)H"),
                                Fragment([Atom("O", "C", 3, 4, 5),
                                          Atom("H", "D", 6, 7, 8),
                                          Atom("H", "D", 9, 1, 2)
                                          ], "water", 0, 1, "O(H)H")
                                ]), molecules)

        self.test_passed = True

    def test_no_settings(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        molecules = molecule_parser.xyz_to_molecules(TestMoleculeParser.monomer_path)

        self.assertEqual(len(molecules), 2)

        for molecule in molecules:
            self.assertEqual(molecule.get_charge(), 0)
            self.assertEqual(molecule.get_spin_multiplicity(), 1)
            self.assertEqual(molecule.get_num_fragments(), 1)
            self.assertEqual(molecule.get_num_atoms(), 3)
            self.assertEqual(molecule.get_symmetry(), "A1B1C1")

        self.assertEqual(molecules[0].get_atoms(), [Atom("O", "A", 1, 2, 3),
                                                    Atom("H", "B", 4, 5, 6),
                                                    Atom("H", "C", 7, 8, 9)
                                                    ])

        self.assertEqual(molecules[1].get_atoms(), [Atom("O", "A", 2, 3, 4),
                                                    Atom("H", "B", 5, 6, 7),
                                                    Atom("H", "C", 8, 9, 10)
                                                    ])

        self.test_passed = True

    def test_confirm_identical_results(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        molecules1 = molecule_parser.xyz_to_molecules(TestMoleculeParser.monomer_path, settings=TestMoleculeParser.monomer_settings)
        molecules2 = list(molecule_parser.parse_training_set_file(TestMoleculeParser.monomer_path, settings=TestMoleculeParser.monomer_settings))

        self.assertEqual(molecules1, molecules2)

        molecules1 = molecule_parser.xyz_to_molecules(TestMoleculeParser.dimer_path, settings=TestMoleculeParser.dimer_settings)
        molecules2 = list(molecule_parser.parse_training_set_file(TestMoleculeParser.dimer_path, settings=TestMoleculeParser.dimer_settings))

        self.assertEqual(molecules1, molecules2)

        molecules1 = molecule_parser.xyz_to_molecules(TestMoleculeParser.trimer_path, settings=TestMoleculeParser.trimer_settings)
        molecules2 = list(molecule_parser.parse_training_set_file(TestMoleculeParser.trimer_path, settings=TestMoleculeParser.trimer_settings))

        self.assertEqual(molecules1, molecules2)

        molecules1 = molecule_parser.xyz_to_molecules(TestMoleculeParser.trimer_path)
        molecules2 = list(molecule_parser.parse_training_set_file(TestMoleculeParser.trimer_path))

        self.assertEqual(molecules1, molecules2)

        self.test_passed = True


suite = unittest.TestLoader().loadTestsFromTestCase(TestMoleculeParser)
