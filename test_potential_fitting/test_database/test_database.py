import unittest, os, psycopg2, random

from potential_fitting.database import Database
from potential_fitting.exceptions import InvalidValueError
from potential_fitting.molecule import Atom, Fragment, Molecule


class TestDatabase(unittest.TestCase):

    @staticmethod
    def get_water_monomer():
        H1 = Atom("H", "A", random.random(), random.random(), random.random())
        H2 = Atom("H", "A", random.random(), random.random(), random.random())
        O1 = Atom("O", "B", random.random(), random.random(), random.random())

        frag1 = Fragment([H1, H2, O1], "H2O", 0, 1, "H1.HO1")

        molecule = Molecule([frag1])

        return molecule

    @staticmethod
    def get_I_monomer():
        I = Atom("I", "A", random.random(), random.random(), random.random())

        frag1 = Fragment([I], "I", 0, 1, "I")

        molecule = Molecule([frag1])

        return molecule

    @staticmethod
    def get_water_dimer():
        H1 = Atom("H", "A", random.random(), random.random(), random.random())
        H2 = Atom("H", "A", random.random(), random.random(), random.random())
        O1 = Atom("O", "B", random.random(), random.random(), random.random())

        frag1 = Fragment([H1, H2, O1], "H2O", 0, 1, "H1.HO1")

        H1 = Atom("H", "A", random.random(), random.random(), random.random())
        H2 = Atom("H", "A", random.random(), random.random(), random.random())
        O1 = Atom("O", "B", random.random(), random.random(), random.random())

        frag2 = Fragment([H1, H2, O1], "H2O", 0, 1, "H1.HO1")

        molecule = Molecule([frag1, frag2])

        return molecule

    @staticmethod
    def get_I_water_dimer():

        I = Atom("I", "A", random.random(), random.random(), random.random())

        frag1 = Fragment([I], "I", 0, 1, "I")

        H1 = Atom("H", "B", random.random(), random.random(), random.random())
        H2 = Atom("H", "B", random.random(), random.random(), random.random())
        O1 = Atom("O", "C", random.random(), random.random(), random.random())

        frag2 = Fragment([H1, H2, O1], "H2O", 0, 1, "H1.HO1")

        molecule = Molecule([frag1, frag2])

        return molecule


    def setUp(self):
        self.config = os.path.join(os.path.dirname(os.path.abspath(__file__)), "local.ini")
        self.database = Database(self.config)
        self.database.annihilate(confirm="confirm")

    def tearDown(self):
        self.database.close()

    def test_set_and_get_batch_size(self):
        self.database.set_batch_size(8)
        self.assertEqual(self.database.get_batch_size(), 8)

        self.database.set_batch_size(15)
        self.assertEqual(self.database.get_batch_size(), 15)

        with self.assertRaises(InvalidValueError):
            self.database.set_batch_size(-1)

        self.assertEqual(self.database.get_batch_size(), 15)

    def test_create(self):
        pass

    def test_execute(self):

        with self.assertRaises(psycopg2.InternalError):
            self.database.execute("RAISE EXCEPTION %s;", ("EXECUTE WORKED",))

    def test_create_postgres_array(self):
        self.assertEqual(self.database.create_postgres_array("A", "B", "C", "D"), "{A,B,C,D}")
        self.assertEqual(self.database.create_postgres_array(1, 2, 3, 4), "{1,2,3,4}")

    def test_add_calculation_and_get_all_calculations(self):

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 0))
        self.assertEqual(len(calculations), 0)

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 10))
        self.assertEqual(len(calculations), 0)

        molecules = []
        for i in range(100):
            molecules.append(self.get_water_monomer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "database_test")

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 100))
        self.assertEqual(len(calculations), 100)

        mols = [calc[0] for calc in calculations]

        molecules = [molecule.get_standard_copy() for molecule in molecules]

        for molecule in molecules:
            self.assertIn((molecule, "testmethod", "testbasis", True, False, [0]), calculations)

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 100))
        self.assertEqual(len(calculations), 0)

        molecules = []
        for i in range(100):
            molecules.append(self.get_water_monomer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 100))
        self.assertEqual(len(calculations), 100)

        molecules = [molecule.get_standard_copy() for molecule in molecules]

        for molecule in molecules:
            self.assertIn((molecule, "testmethod", "testbasis", False, False, [0]), calculations)

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 100))
        self.assertEqual(len(calculations), 0)

        molecules = []
        for i in range(100):
            molecules.append(self.get_water_dimer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "database_test")

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 500))
        self.assertEqual(len(calculations), 500)

        molecules = [molecule.get_standard_copy() for molecule in molecules]

        for molecule in molecules:
            self.assertIn((molecule, "testmethod", "testbasis", True, False, [0]), calculations)
            self.assertIn((molecule, "testmethod", "testbasis", True, False, [1]), calculations)
            self.assertIn((molecule, "testmethod", "testbasis", True, True, [0]), calculations)
            self.assertIn((molecule, "testmethod", "testbasis", True, True, [1]), calculations)
            self.assertIn((molecule, "testmethod", "testbasis", True, False, [0, 1]), calculations)

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 100))
        self.assertEqual(len(calculations), 0)

        molecules = []
        for i in range(100):
            molecules.append(self.get_water_dimer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 500))
        self.assertEqual(len(calculations), 300)

        molecules = [molecule.get_standard_copy() for molecule in molecules]

        for molecule in molecules:
            self.assertIn((molecule, "testmethod", "testbasis", False, False, [0]), calculations)
            self.assertIn((molecule, "testmethod", "testbasis", False, False, [1]), calculations)
            self.assertIn((molecule, "testmethod", "testbasis", False, False, [0, 1]), calculations)

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 200))
        self.assertEqual(len(calculations), 0)

    def test_delete_calculations(self):

        molecules = []
        for i in range(100):
            molecules.append(self.get_water_monomer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "database_test")

        self.database.delete_calculations(molecules, "testmethod", "testbasis", True, "database_test", delete_complete_calculations = False)

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 100))
        self.assertEqual(len(calculations), 0)

        opt_mol = self.get_water_monomer()
        opt_energy = random.random()
        molecules = []
        for i in range(100):
            molecules.append(self.get_water_monomer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", False, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do=101)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                energy = random.random()

            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, False, energy, "some log test"])

        self.database.set_properties(calculation_results)

        self.database.delete_calculations(molecules, "testmethod", "testbasis", False, "database_test", delete_complete_calculations = False)
        self.database.delete_calculations([opt_mol], "testmethod", "testbasis", False, "database_test", delete_complete_calculations = False)

        with self.assertRaises(psycopg2.InternalError):
            training_set = list(self.database.get_1B_training_set("H2O", ["H2O"], ["H1.HO1"], "testmethod", "testbasis", True, "database_test"))

    def test_delete_all_calculations(self):

        molecules = []
        for i in range(100):
            molecules.append(self.get_water_monomer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "database_test")

        self.database.delete_all_calculations("H2O", "testmethod", "testbasis", True, "database_test", delete_complete_calculations = False)

        calculations = list(self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 100))
        self.assertEqual(len(calculations), 0)

        opt_mol = self.get_water_monomer()
        opt_energy = random.random()
        molecules = []
        for i in range(100):
            molecules.append(self.get_water_monomer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", False, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do=101)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                energy = random.random()

            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, False, energy, "some log test"])

        self.database.set_properties(calculation_results)

        self.database.delete_all_calculations("H2O", "testmethod", "testbasis", False, "database_test", delete_complete_calculations = False)

        with self.assertRaises(psycopg2.InternalError):
            training_set = list(self.database.get_1B_training_set("H2O", ["H2O"], ["H1.HO1"], "testmethod", "testbasis", True, "database_test"))

    def test_set_properties_and_get_1B_training_set(self):
        opt_mol = self.get_water_monomer()
        opt_energy = random.random()
        molecules = []
        for i in range(100):
            molecules.append(self.get_water_monomer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", True, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 101)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                energy = random.random()

            calculation_results.append([molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(self.database.get_1B_training_set("H2O", ["H2O"], ["H1.HO1"], "testmethod", "testbasis", True, "database_test"))
        self.assertEqual(len(training_set), 101)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index] = tuple(training_set[index])

        for molecule, method, basis, cp, use_cp, frag_indices, result, energy, log_text in calculation_results:
            self.assertIn((molecule.get_reorder_copy(["H2O"], ["H1.HO1"]), round(energy - opt_energy, 5)), training_set)

    def test_set_properties_and_get_2B_training_set(self):

        # water_water dimer

        opt_mol = self.get_water_monomer()
        opt_energy = random.random()

        molecules = []
        energies = []
        for i in range(100):
            molecules.append(self.get_water_dimer())
            energies.append([random.random(), random.random(), random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", False, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 301)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                index = molecules.index(molecule.get_reorder_copy(["H2O", "H2O"], ["H1.HO1", "H1.HO1"]))
                if frag_indices == [0]:
                    energy = energies[index][0]
                elif frag_indices == [1]:
                    energy = energies[index][1]
                else:
                    energy = energies[index][2]
            calculation_results.append([molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(self.database.get_2B_training_set("H2O-H2O", ["H2O", "H2O"], ["H1.HO1", "H1.HO1"], "testmethod", "testbasis", False, "database_test"))

        self.assertEqual(len(training_set), 100)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3] = round(training_set[index][3], 5)
            training_set[index][4] = round(training_set[index][4], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            monomer1 = energies[index][0] - opt_energy
            monomer2 = energies[index][1] - opt_energy
            interaction = energies[index][2] - energies[index][1] - energies[index][0]
            binding = monomer1 + monomer2 + interaction
            self.assertIn((molecule, round(binding, 5), round(interaction, 5), round(monomer1, 5), round(monomer2, 5)), training_set)

        # I_water dimer

        opt_I = self.get_I_monomer()
        opt_I_energy = random.random()

        opt_water = self.get_water_monomer()
        opt_water_energy = random.random()

        molecules = []
        energies = []
        for i in range(10):
            molecules.append(self.get_I_water_dimer())
            energies.append([random.random(), random.random(), random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test_2")
        self.database.add_calculations([opt_I], "testmethod", "testbasis", False, "database_test_2", optimized=True)
        self.database.add_calculations([opt_water], "testmethod", "testbasis", False, "database_test_2", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test_2", calculations_to_do = 301)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_I.get_standard_copy():
                energy = opt_I_energy
            elif molecule == opt_water.get_standard_copy():
                energy = opt_water_energy
            else:
                index = molecules.index(molecule.get_reorder_copy(["I", "H2O"], ["I", "H1.HO1"]))
                if frag_indices == [0]:
                    energy = energies[index][1]
                elif frag_indices == [1]:
                    energy = energies[index][0]
                else:
                    energy = energies[index][2]
            calculation_results.append([molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(self.database.get_2B_training_set("H2O-I", ["I", "H2O"], ["I", "H1.HO1"], "testmethod", "testbasis", False, "database_test_2"))
        self.assertEqual(len(training_set), 10)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3] = round(training_set[index][3], 5)
            training_set[index][4] = round(training_set[index][4], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            monomer1 = energies[index][0] - opt_I_energy
            monomer2 = energies[index][1] - opt_water_energy
            interaction = energies[index][2] - energies[index][1] - energies[index][0]
            binding = monomer1 + monomer2 + interaction
            self.assertIn((molecule, round(binding, 5), round(interaction, 5), round(monomer1, 5), round(monomer2, 5)), training_set)

    def test_import_calculations(self):

        # 1B

        opt_mol = self.get_water_monomer()
        opt_energy = random.random()
        molecule_energies_pairs = []
        for i in range(100):
            molecule_energies_pairs.append((self.get_water_monomer(), [random.random()]))

        self.database.import_calculations(molecule_energies_pairs, "testmethod", "testbasis", False, "database_test", optimized=False)
        self.database.import_calculations([(opt_mol, [opt_energy])], "testmethod", "testbasis", False, "database_test", optimized=True)

        training_set = list(self.database.get_1B_training_set("H2O", ["H2O"], ["H1.HO1"], "testmethod", "testbasis", False, "database_test"))
        self.assertEqual(len(training_set), 101)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index] = tuple(training_set[index])

        for molecule, energies in molecule_energies_pairs:
            self.assertIn((molecule.get_reorder_copy(["H2O"], ["H1.HO1"]), round(energies[0] - opt_energy, 5)), training_set)

        # 2B water dimer

        molecule_energies_pairs = []
        for i in range(100):
            molecule_energies_pairs.append((self.get_water_dimer(), [random.random(), random.random(), random.random()]))

        self.database.import_calculations(molecule_energies_pairs, "testmethod", "testbasis", False, "database_test", optimized=False)
        self.database.import_calculations([(opt_mol, [opt_energy])], "testmethod", "testbasis", False, "database_test", optimized=True)

        training_set = list(self.database.get_2B_training_set("H2O-H2O", ["H2O", "H2O"], ["H1.HO1", "H1.HO1"], "testmethod", "testbasis", False, "database_test"))
        self.assertEqual(len(training_set), 100)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3] = round(training_set[index][3], 5)
            training_set[index][4] = round(training_set[index][4], 5)
            training_set[index] = tuple(training_set[index])

        for molecule, energies in molecule_energies_pairs:
            monomer1 = energies[0] - opt_energy
            monomer2 = energies[1] - opt_energy
            interaction = energies[2] - energies[1] - energies[0]
            binding = monomer1 + monomer2 + interaction
            self.assertIn((molecule, round(binding, 5), round(interaction, 5), round(monomer1, 5), round(monomer2, 5)),
                          training_set)

        # 2B I_water dimer

        opt_I = self.get_I_monomer()
        opt_I_energy = random.random()

        molecule_energies_pairs = []
        for i in range(10):
            molecule_energies_pairs.append(
                (self.get_I_water_dimer(), [random.random(), random.random(), random.random()]))

        self.database.import_calculations(molecule_energies_pairs, "testmethod", "testbasis", False, "database_test",
                                          optimized=False)
        self.database.import_calculations([(opt_mol, [opt_energy])], "testmethod", "testbasis", False, "database_test",
                                          optimized=True)
        self.database.import_calculations([(opt_I, [opt_I_energy])], "testmethod", "testbasis", False, "database_test",
                                          optimized=True)

        training_set = list(
            self.database.get_2B_training_set("H2O-I", ["I", "H2O"], ["I", "H1.HO1"], "testmethod",
                                              "testbasis", False, "database_test"))
        self.assertEqual(len(training_set), 10)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3] = round(training_set[index][3], 5)
            training_set[index][4] = round(training_set[index][4], 5)
            training_set[index] = tuple(training_set[index])

        for molecule, energies in molecule_energies_pairs:
            monomer1 = energies[0] - opt_I_energy
            monomer2 = energies[1] - opt_energy
            interaction = energies[2] - energies[1] - energies[0]
            binding = monomer1 + monomer2 + interaction
            self.assertIn((molecule, round(binding, 5), round(interaction, 5), round(monomer1, 5), round(monomer2, 5)),
                          training_set)

    def test_export_calculations(self):
        # 1B

        molecule_energies_pairs = []
        for i in range(10):
            molecule_energies_pairs.append((self.get_water_monomer(), [random.random()]))

        self.database.import_calculations(molecule_energies_pairs, "testmethod", "testbasis", False, "database_test",
                                          optimized=False)

        calculations = list(
            self.database.export_calculations(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                              "database_test"))
        self.assertEqual(len(calculations), 10)

        for index in range(len(calculations)):
            calculations[index] = list(calculations[index])
            calculations[index][1][0] = round(calculations[index][1][0], 5)
            calculations[index] = tuple(calculations[index])

        for molecule, energies in molecule_energies_pairs:
            self.assertIn((molecule.get_reorder_copy(["H2O"], ["H1.HO1"]), [round(energies[0], 5)]),
            calculations)

        # 2B
        molecule_energies_pairs = []
        for i in range(10):
            molecule_energies_pairs.append((self.get_I_water_dimer(), [random.random(), random.random(), random.random(), random.random(), random.random()]))


        self.database.import_calculations(molecule_energies_pairs, "testmethod", "testbasis", True,
                                          "database_test",
                                          optimized=False)

        calculations = list(
            self.database.export_calculations(["I", "H2O"], ["I", "H1.HO1"], "testmethod", "testbasis", True,
                                              "database_test"))
        self.assertEqual(len(calculations), 10)

        for index in range(len(calculations)):
            calculations[index] = list(calculations[index])
            calculations[index][1][0] = round(calculations[index][1][0], 5)
            calculations[index][1][1] = round(calculations[index][1][1], 5)
            calculations[index][1][2] = round(calculations[index][1][2], 5)
            calculations[index][1][3] = round(calculations[index][1][3], 5)
            calculations[index][1][4] = round(calculations[index][1][4], 5)
            calculations[index] = tuple(calculations[index])

        for molecule, energies in molecule_energies_pairs:
            self.assertIn((molecule.get_reorder_copy(["I", "H2O"], ["I", "H1.HO1"]), [round(energies[0], 5), round(energies[1], 5), round(energies[2], 5), round(energies[3], 5), round(energies[4], 5)]),
                          calculations)


suite = unittest.TestLoader().loadTestsFromTestCase(TestDatabase)