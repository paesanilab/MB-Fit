import unittest, os, random

from test_mbfit.test_case_with_id import TestCaseWithId
from mbfit.database import Database
from mbfit.exceptions import InvalidValueError, DatabaseConnectionError, DatabaseOperationError
from mbfit.molecule import Atom, Fragment, Molecule

# only import psycopg2 if it is installed.
try:
	import psycopg2
except ModuleNotFoundError:
    pass

def psycopg2_installed():
    try:
        import psycopg2
        return True
    except ModuleNotFoundError:
        return False

def local_db_installed():
    try:
        config = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_user1.ini")
        database = Database(config)
        database.close()
        return True
    except DatabaseConnectionError:
        return False

@unittest.skipUnless(psycopg2_installed() and local_db_installed(),"psycopg2 or a local test database is not installed, so database cannot be tested.")
class TestDatabase(TestCaseWithId):
    def __init__(self, *args, **kwargs):
        super(TestDatabase, self).__init__(*args, **kwargs)
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

    @staticmethod
    def get_water_monomer():
        H1 = Atom("H", "A", random.random(), random.random(), random.random())
        H2 = Atom("H", "A", random.random(), random.random(), random.random())
        O1 = Atom("O", "B", random.random(), random.random(), random.random())

        frag1 = Fragment([H1, H2, O1], "H2O", 0, 1, "H1.HO1")

        molecule = Molecule([frag1])

        return molecule

    @staticmethod
    def get_ethane_monomer():

        frag1 = Fragment([Atom("C", "A", random.random(), random.random(), random.random()),
                          Atom("H", "B", random.random(), random.random(), random.random()),
                          Atom("H", "B", random.random(), random.random(), random.random()),
                          Atom("H", "B", random.random(), random.random(), random.random()),
                          Atom("C", "A", random.random(), random.random(), random.random()),
                          Atom("H", "B", random.random(), random.random(), random.random()),
                          Atom("H", "B", random.random(), random.random(), random.random()),
                          Atom("H", "B", random.random(), random.random(), random.random())], "ethane", 0, 1, "C1(H)(H)H.C1(H)(H)H")

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
    def get_water_trimer():
        H1 = Atom("H", "A", random.random(), random.random(), random.random())
        H2 = Atom("H", "A", random.random(), random.random(), random.random())
        O1 = Atom("O", "B", random.random(), random.random(), random.random())

        frag1 = Fragment([H1, H2, O1], "H2O", 0, 1, "H1.HO1")

        H1 = Atom("H", "A", random.random(), random.random(), random.random())
        H2 = Atom("H", "A", random.random(), random.random(), random.random())
        O1 = Atom("O", "B", random.random(), random.random(), random.random())

        frag2 = Fragment([H1, H2, O1], "H2O", 0, 1, "H1.HO1")

        H1 = Atom("H", "A", random.random(), random.random(), random.random())
        H2 = Atom("H", "A", random.random(), random.random(), random.random())
        O1 = Atom("O", "B", random.random(), random.random(), random.random())

        frag3 = Fragment([H1, H2, O1], "H2O", 0, 1, "H1.HO1")

        molecule = Molecule([frag1, frag2, frag3])

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
        self.config = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_user1.ini")
        self.config2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_user2.ini")
        self.config3 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_user3.ini")

        self.database = Database(self.config)
        self.database2 = Database(self.config2)
        self.database3 = Database(self.config3)

        self.database.annihilate(confirm="confirm")
        self.database.save()

    def tearDown(self):
        self.database2.close()
        self.database3.close()
        self.database.annihilate(confirm="confirm")
        self.database.save()
        self.database.close()

        super().tearDown()

    def test_set_and_get_batch_size(self):

        self.database.set_batch_size(8)
        self.assertEqual(self.database.get_batch_size(), 8)

        self.database.set_batch_size(15)
        self.assertEqual(self.database.get_batch_size(), 15)

        with self.assertRaises(InvalidValueError):
            self.database.set_batch_size(-1)

        self.assertEqual(self.database.get_batch_size(), 15)

        self.test_passed = True

    def test_create(self):

        self.test_passed = True

    def test_execute(self):

        with self.assertRaises(DatabaseOperationError):
            self.database.execute("RAISE EXCEPTION %s;", ("EXECUTE WORKED",))

        self.test_passed = True

    def test_create_postgres_array(self):
        self.assertEqual(self.database.create_postgres_array("A", "B", "C", "D"), "{A,B,C,D}")
        self.assertEqual(self.database.create_postgres_array(1, 2, 3, 4), "{1,2,3,4}")

        self.test_passed = True

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

        self.test_passed = True

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

        with self.assertRaises(DatabaseOperationError):
            training_set = list(self.database.get_1B_training_set("H2O", ["H2O"], ["H1.HO1"], "testmethod", "testbasis", True, "database_test"))

        self.test_passed = True

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

        with self.assertRaises(DatabaseOperationError):
            training_set = list(self.database.get_1B_training_set("H2O", ["H2O"], ["H1.HO1"], "testmethod", "testbasis", True, "database_test"))

        self.test_passed = True

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

        self.test_passed = True

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
            nb_energy = energies[index][2] - energies[index][1] - energies[index][0]
            binding = monomer1 + monomer2 + nb_energy
            self.assertIn((molecule, round(binding, 5), round(nb_energy, 5), round(monomer1, 5), round(monomer2, 5)), training_set)

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
        self.database.add_calculations([opt_water], "testmethod", "testbasis", False, "database_test_2",
                                       optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test_2", calculations_to_do=301)

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
            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(
            self.database.get_2B_training_set("H2O-I", ["I", "H2O"], ["I", "H1.HO1"], "testmethod", "testbasis",
                                              False, "database_test_2"))
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
            nb_energy = energies[index][2] - energies[index][1] - energies[index][0]
            binding = monomer1 + monomer2 + nb_energy
            self.assertIn(
                (molecule, round(binding, 5), round(nb_energy, 5), round(monomer1, 5), round(monomer2, 5)),
                training_set)

        self.test_passed = True

    def test_set_properties_and_get_training_set_1B(self):

        # no cp

        opt_mol = self.get_water_monomer()
        opt_energy = random.random()

        molecules = []
        energies = []
        for i in range(100):
            molecules.append(self.get_water_monomer())
            energies.append([random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", False, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do=101)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                index = molecules.index(molecule.get_reorder_copy(["H2O"], ["H1.HO1"]))
                if frag_indices == [0]:
                    energy = energies[index][0]
                elif frag_indices == [1]:
                    energy = energies[index][1]
                else:
                    energy = energies[index][2]
            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(
            self.database.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))

        self.assertEqual(len(training_set), 101)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3][0] = round(training_set[index][3][0], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            binding = energies[index][0] - opt_energy
            self.assertIn((molecule, round(binding, 5), round(binding, 5), [round(binding, 5)]), training_set)

        self.test_passed = True

    def test_set_properties_and_get_training_set_2B(self):

        # no cp

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

        training_set = list(self.database.get_training_set(["H2O", "H2O"], ["H1.HO1", "H1.HO1"], "testmethod", "testbasis", False, "database_test"))

        self.assertEqual(len(training_set), 100)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3][0] = round(training_set[index][3][0], 5)
            training_set[index][3][1] = round(training_set[index][3][1], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            monomer1 = energies[index][0] - opt_energy
            monomer2 = energies[index][1] - opt_energy
            nb_energy = energies[index][2] - energies[index][1] - energies[index][0]
            binding = monomer1 + monomer2 + nb_energy
            self.assertIn((molecule, round(binding, 5), round(nb_energy, 5), [round(monomer1, 5), round(monomer2, 5)]), training_set)

        # yes cp

        molecules = []
        energies = []
        for i in range(100):
            molecules.append(self.get_water_dimer())
            energies.append([random.random(), random.random(), random.random(), random.random(), random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", True, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do=501)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                index = molecules.index(molecule.get_reorder_copy(["H2O", "H2O"], ["H1.HO1", "H1.HO1"]))
                if frag_indices == [0]:
                    if use_cp:
                        energy = energies[index][0]
                    else:
                        energy = energies[index][1]
                elif frag_indices == [1]:
                    if use_cp:
                        energy = energies[index][2]
                    else:
                        energy = energies[index][3]
                else:
                    energy = energies[index][4]
            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(
            self.database.get_training_set(["H2O", "H2O"], ["H1.HO1", "H1.HO1"], "testmethod", "testbasis", True,
                                           "database_test"))

        self.assertEqual(len(training_set), 100)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3][0] = round(training_set[index][3][0], 5)
            training_set[index][3][1] = round(training_set[index][3][1], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            monomer1 = energies[index][1] - opt_energy
            monomer2 = energies[index][3] - opt_energy
            nb_energy = energies[index][4] - energies[index][0] - energies[index][2]
            binding = monomer1 + monomer2 + nb_energy

            self.assertIn((molecule, round(binding, 5), round(nb_energy, 5), [round(monomer1, 5), round(monomer2, 5)]), training_set)

        self.test_passed = True

    def test_set_properties_and_get_training_set_3B(self):

        # no cp

        opt_mol = self.get_water_monomer()
        opt_energy = random.random()

        molecules = []
        energies = []
        for i in range(100):
            molecules.append(self.get_water_trimer())
            energies.append([random.random(), random.random(), random.random(), random.random(), random.random(), random.random(), random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", False, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do = 701)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                index = molecules.index(molecule.get_reorder_copy(["H2O", "H2O", "H2O"], ["H1.HO1", "H1.HO1", "H1.HO1"]))
                if frag_indices == [0]:
                    energy = energies[index][0]
                elif frag_indices == [1]:
                    energy = energies[index][1]
                elif frag_indices == [2]:
                    energy = energies[index][2]
                elif frag_indices == [0, 1]:
                    energy = energies[index][3]
                elif frag_indices == [0, 2]:
                    energy = energies[index][4]
                elif frag_indices == [1, 2]:
                    energy = energies[index][5]
                else:
                    energy = energies[index][6]
            calculation_results.append([molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(self.database.get_training_set(["H2O", "H2O", "H2O"], ["H1.HO1", "H1.HO1", "H1.HO1"], "testmethod", "testbasis", False, "database_test"))

        self.assertEqual(len(training_set), 100)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3][0] = round(training_set[index][3][0], 5)
            training_set[index][3][1] = round(training_set[index][3][1], 5)
            training_set[index][3][2] = round(training_set[index][3][2], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            monomer1 = energies[index][0] - opt_energy
            monomer2 = energies[index][1] - opt_energy
            monomer3 = energies[index][2] - opt_energy
            nb_energy = energies[index][6] - energies[index][5] - energies[index][4] - energies[index][3] + energies[index][2] + energies[index][1] + energies[index][0]
            binding = monomer1 + monomer2 + monomer3 + nb_energy
            self.assertIn((molecule, round(binding, 5), round(nb_energy, 5), [round(monomer1, 5), round(monomer2, 5), round(monomer3, 5)]), training_set)

        # yes cp

        molecules = []
        energies = []
        for i in range(100):
            molecules.append(self.get_water_trimer())
            energies.append([random.random(), random.random(), random.random(), random.random(),
                             random.random(), random.random(), random.random(), random.random(),
                             random.random(), random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", True, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do=1301)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                index = molecules.index(
                    molecule.get_reorder_copy(["H2O", "H2O", "H2O"], ["H1.HO1", "H1.HO1", "H1.HO1"]))
                if frag_indices == [0]:
                    if use_cp:
                        energy = energies[index][0]
                    else:
                        energy = energies[index][1]
                elif frag_indices == [1]:
                    if use_cp:
                        energy = energies[index][2]
                    else:
                        energy = energies[index][3]
                elif frag_indices == [2]:
                    if use_cp:
                        energy = energies[index][4]
                    else:
                        energy = energies[index][5]
                elif frag_indices == [0, 1]:
                    energy = energies[index][6]
                elif frag_indices == [0, 2]:
                    energy = energies[index][7]
                elif frag_indices == [1, 2]:
                    energy = energies[index][8]
                else:
                    energy = energies[index][9]
            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(
            self.database.get_training_set(["H2O", "H2O", "H2O"], ["H1.HO1", "H1.HO1", "H1.HO1"], "testmethod",
                                           "testbasis", True, "database_test"))

        self.assertEqual(len(training_set), 100)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3][0] = round(training_set[index][3][0], 5)
            training_set[index][3][1] = round(training_set[index][3][1], 5)
            training_set[index][3][2] = round(training_set[index][3][2], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            monomer1 = energies[index][1] - opt_energy
            monomer2 = energies[index][3] - opt_energy
            monomer3 = energies[index][5] - opt_energy
            nb_energy = energies[index][9] - energies[index][8] - energies[index][7] - energies[index][6] + \
                          energies[index][4] + energies[index][2] + energies[index][0]
            binding = monomer1 + monomer2 + monomer3 + nb_energy
            self.assertIn((molecule, round(binding, 5), round(nb_energy, 5), [round(monomer1, 5), round(monomer2, 5), round(monomer3, 5)]), training_set)

        self.test_passed = True

    def test_set_properties_and_get_training_set_nested_symmetry(self):

        # no cp

        opt_mol = self.get_ethane_monomer()
        opt_energy = random.random()

        molecules = []
        energies = []
        for i in range(5):
            molecules.append(self.get_ethane_monomer())
            energies.append([random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", False, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do=6)

        calculation_results = []
        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                index = molecules.index(molecule.get_reorder_copy(["ethane"], ["C1(H)(H)H.C1(H)(H)H"]))
                if frag_indices == [0]:
                    energy = energies[index][0]
                elif frag_indices == [1]:
                    energy = energies[index][1]
                else:
                    energy = energies[index][2]
            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(
            self.database.get_training_set(["ethane"], ["C1(H)(H)H.C1(H)(H)H"], "testmethod", "testbasis", False,
                                           "database_test"))

        self.assertEqual(len(training_set), 6)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3][0] = round(training_set[index][3][0], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            binding = energies[index][0] - opt_energy
            self.assertIn((molecule, round(binding, 5), round(binding, 5), [round(binding, 5)]), training_set)

        self.test_passed = True

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
            nb_energy = energies[2] - energies[1] - energies[0]
            binding = monomer1 + monomer2 + nb_energy
            self.assertIn((molecule, round(binding, 5), round(nb_energy, 5), round(monomer1, 5), round(monomer2, 5)),
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
            nb_energy = energies[2] - energies[1] - energies[0]
            binding = monomer1 + monomer2 + nb_energy
            self.assertIn((molecule, round(binding, 5), round(nb_energy, 5), round(monomer1, 5), round(monomer2, 5)),
                          training_set)

        self.test_passed = True

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

        self.test_passed = True

    def test_read_privileges(self):

        # First, create the training set owned by test_user1

        opt_mol = self.get_water_monomer()
        opt_energy = random.random()

        molecules = []
        energies = []
        for i in range(100):
            molecules.append(self.get_water_monomer())
            energies.append([random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", False, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do=101)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                index = molecules.index(molecule.get_reorder_copy(["H2O"], ["H1.HO1"]))
                if frag_indices == [0]:
                    energy = energies[index][0]
                elif frag_indices == [1]:
                    energy = energies[index][1]
                else:
                    energy = energies[index][2]
            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(
            self.database.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))

        self.assertEqual(len(training_set), 101)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3][0] = round(training_set[index][3][0], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            binding = energies[index][0] - opt_energy
            self.assertIn((molecule, round(binding, 5), round(binding, 5), [round(binding, 5)]), training_set)

        self.database.save()

        # Now, test_user2 tries to read the training set without permissions.

        with self.assertRaises(DatabaseOperationError):
            try:
                list(self.database2.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have read privileges on training set database_test", str(e))
                raise

        with self.assertRaises(DatabaseOperationError):
            try:
                list(self.database2.export_calculations(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have read privileges on training set database_test", str(e))
                raise

        with self.assertRaises(DatabaseOperationError):
            try:
                list(self.database2.get_failed("H2O", ["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have read privileges on training set database_test", str(e))
                raise

        # test_user2 is given read privileges, and can now perform read operations.
        self.database.grant_read_privilege("test_user2", "database_test")
        self.database.save()

        training_set = list(self.database2.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                   "database_test"))
        self.assertEqual(len(training_set), 101)

        calculations = list(self.database2.export_calculations(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))
        self.assertEqual(len(calculations), 101)

        failed = list(self.database2.get_failed("H2O", ["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))
        self.assertEqual(len(failed), 0)

        # test_user2 gets their read privileges revoked, and can no longer perform read operations.

        self.database.revoke_read_privilege("test_user2", "database_test")
        self.database.save()

        with self.assertRaises(DatabaseOperationError):
            try:
                list(
                    self.database2.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                                    "database_test"))
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have read privileges on training set database_test", str(e))
                raise

        with self.assertRaises(DatabaseOperationError):
            try:
                list(self.database2.export_calculations(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                                        "database_test"))
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have read privileges on training set database_test", str(e))
                raise

        with self.assertRaises(DatabaseOperationError):
            try:
                list(self.database2.get_failed("H2O", ["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have read privileges on training set database_test", str(e))
                raise

        self.test_passed = True

    def test_write_privileges(self):

        # First, create the training set owned by test_user1

        opt_mol = self.get_water_monomer()
        opt_energy = random.random()

        molecules = []
        energies = []
        for i in range(100):
            molecules.append(self.get_water_monomer())
            energies.append([random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", False, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do=101)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                index = molecules.index(molecule.get_reorder_copy(["H2O"], ["H1.HO1"]))
                if frag_indices == [0]:
                    energy = energies[index][0]
                elif frag_indices == [1]:
                    energy = energies[index][1]
                else:
                    energy = energies[index][2]
            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(
            self.database.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))

        self.assertEqual(len(training_set), 101)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3][0] = round(training_set[index][3][0], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            binding = energies[index][0] - opt_energy
            self.assertIn((molecule, round(binding, 5), round(binding, 5), [round(binding, 5)]), training_set)

        self.database.save()

        # Now, test_user2 tries to write to the training set without permissions.

        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have write privileges on training set database_test", str(e))
                raise

        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.import_calculations([(molecule, energy) for molecule, energy in zip(molecules, energies)], "testmethod", "testbasis", False, "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have write privileges on training set database_test", str(e))
                raise

        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.delete_calculations(molecules, "testmethod", "testbasis", False, "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have write privileges on training set database_test", str(e))
                raise

        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.delete_all_calculations("H2O", "testmethod", "testbasis", False, "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have write privileges on training set database_test", str(e))
                raise

        # test_user2 is given write privileges, and can now perform write operations.
        self.database.grant_write_privilege("test_user2", "database_test")
        self.database.save()

        molecules = []
        energies = []
        for i in range(10):
            molecules.append(self.get_water_monomer())
            energies.append([random.random()])

        self.database2.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")

        calculations = self.database2.get_all_calculations("testclient", "database_test", calculations_to_do=10)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                index = molecules.index(molecule.get_reorder_copy(["H2O"], ["H1.HO1"]))
                if frag_indices == [0]:
                    energy = energies[index][0]
                elif frag_indices == [1]:
                    energy = energies[index][1]
                else:
                    energy = energies[index][2]
            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database2.set_properties(calculation_results)
        self.database2.save()

        training_set = list(
            self.database.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))

        self.assertEqual(len(training_set), 111)

        molecules = []
        energies = []
        for i in range(10):
            molecules.append(self.get_water_monomer())
            energies.append([random.random()])

        self.database2.import_calculations([(molecule, energy) for molecule, energy in zip(molecules, energies)],
                                           "testmethod", "testbasis", False, "database_test")
        self.database2.save()

        training_set = list(
            self.database.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))

        self.assertEqual(len(training_set), 121)

        self.database2.delete_calculations(molecules, "testmethod", "testbasis", False, "database_test", delete_complete_calculations = False)
        self.database2.save()

        training_set = list(
            self.database.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))

        self.assertEqual(len(training_set), 111)

        self.database2.delete_all_calculations("H2O", "testmethod", "testbasis", False, "database_test", delete_complete_calculations = False)
        self.database2.save()

        # with all calculations deleted, there will be no optimized energy for H2O
        with self.assertRaises(DatabaseOperationError):
            training_set = list(
                self.database.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                               "database_test"))

        self.database2.import_calculations([(molecule, energy) for molecule, energy in zip(molecules, energies)],
                                           "testmethod", "testbasis", False, "database_test")
        self.database2.save()

        # test_user2 gets their write privileges revoked, and can no longer perform write operations.
        self.database.revoke_write_privilege("test_user2", "database_test")
        self.database.save()

        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have write privileges on training set database_test", str(e))
                raise

        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.import_calculations([(molecule, energy) for molecule, energy in zip(molecules, energies)], "testmethod", "testbasis", False, "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have write privileges on training set database_test", str(e))
                raise

        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.delete_calculations(molecules, "testmethod", "testbasis", False, "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have write privileges on training set database_test", str(e))
                raise

        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.delete_all_calculations("H2O", "testmethod", "testbasis", False, "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have write privileges on training set database_test", str(e))
                raise

        self.test_passed = True

    def test_admin_privileges(self):

        # First, create the training set owned by test_user1

        opt_mol = self.get_water_monomer()
        opt_energy = random.random()

        molecules = []
        energies = []
        for i in range(100):
            molecules.append(self.get_water_monomer())
            energies.append([random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "database_test")
        self.database.add_calculations([opt_mol], "testmethod", "testbasis", False, "database_test", optimized=True)

        calculations = self.database.get_all_calculations("testclient", "database_test", calculations_to_do=101)

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            if molecule == opt_mol.get_standard_copy():
                energy = opt_energy
            else:
                index = molecules.index(molecule.get_reorder_copy(["H2O"], ["H1.HO1"]))
                if frag_indices == [0]:
                    energy = energies[index][0]
                elif frag_indices == [1]:
                    energy = energies[index][1]
                else:
                    energy = energies[index][2]
            calculation_results.append(
                [molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = list(
            self.database.get_training_set(["H2O"], ["H1.HO1"], "testmethod", "testbasis", False,
                                           "database_test"))

        self.assertEqual(len(training_set), 101)

        for index in range(len(training_set)):
            training_set[index] = list(training_set[index])
            training_set[index][1] = round(training_set[index][1], 5)
            training_set[index][2] = round(training_set[index][2], 5)
            training_set[index][3][0] = round(training_set[index][3][0], 5)
            training_set[index] = tuple(training_set[index])

        for molecule in molecules:
            index = molecules.index(molecule)
            binding = energies[index][0] - opt_energy
            self.assertIn((molecule, round(binding, 5), round(binding, 5), [round(binding, 5)]), training_set)

        self.database.save()

        # Now, test_user2 tries to change permissions to the training set without admin permissions.

        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.grant_read_privilege("test_user2", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.grant_write_privilege("test_user2", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.grant_admin_privilege("test_user2", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.revoke_read_privilege("test_user1", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.revoke_write_privilege("test_user1", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.revoke_admin_privilege("test_user1", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise

        # test_user2 is given admin privileges, and can now perform admin operations.
        self.database.grant_admin_privilege("test_user2", "database_test")
        self.database.save()

        self.database2.revoke_admin_privilege("test_user1", "database_test")
        self.database2.revoke_write_privilege("test_user1", "database_test")
        self.database2.revoke_read_privilege("test_user1", "database_test")

        # Admins cannot revoke the last admin's privileges.
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.revoke_admin_privilege("test_user2", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 is the last admin of this training set and many not be removed", str(e))
                raise

        self.database2.grant_admin_privilege("test_user1", "database_test")
        self.database2.grant_write_privilege("test_user1", "database_test")
        self.database2.grant_read_privilege("test_user1", "database_test")

        # test_user2 gets their write privileges revoked, and can no longer perform write operations.
        self.database.revoke_admin_privilege("test_user2", "database_test")
        self.database.save()

        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.grant_read_privilege("test_user2", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.grant_write_privilege("test_user2", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.grant_admin_privilege("test_user2", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.revoke_read_privilege("test_user1", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.revoke_write_privilege("test_user1", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise
        with self.assertRaises(DatabaseOperationError):
            try:
                self.database2.revoke_admin_privilege("test_user1", "database_test")
            except DatabaseOperationError as e:
                self.assertIn("User test_user2 does not have admin privileges on training set database_test", str(e))
                raise

        self.test_passed = True

suite = unittest.TestLoader().loadTestsFromTestCase(TestDatabase)
