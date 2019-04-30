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

        frag1 = Fragment("H2O", random.randint(0, 10), random.randint(1, 10))
        frag1.add_atom(H1)
        frag1.add_atom(H2)
        frag1.add_atom(O1)

        molecule = Molecule()
        molecule.add_fragment(frag1)

        return molecule

    @staticmethod
    def get_water_dimer():
        H1 = Atom("H", "A", random.random(), random.random(), random.random())
        H2 = Atom("H", "A", random.random(), random.random(), random.random())
        O1 = Atom("O", "B", random.random(), random.random(), random.random())

        frag1 = Fragment("H2O", random.randint(0, 10), random.randint(1, 10))
        frag1.add_atom(H1)
        frag1.add_atom(H2)
        frag1.add_atom(O1)

        H1 = Atom("H", "A", random.random(), random.random(), random.random())
        H2 = Atom("H", "A", random.random(), random.random(), random.random())
        O1 = Atom("O", "B", random.random(), random.random(), random.random())

        frag2 = Fragment("H2O", random.randint(0, 10), random.randint(1, 10))
        frag2.add_atom(H1)
        frag2.add_atom(H2)
        frag2.add_atom(O1)

        molecule = Molecule()
        molecule.add_fragment(frag1)
        molecule.add_fragment(frag2)

        return molecule

    def setUp(self):
        self.config = os.path.join(os.path.dirname(os.path.abspath(__file__)), "database.ini")
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
        self.database.save()

        pass

    def test_create_postgres_array(self):
        self.assertEqual(self.database.create_postgres_array("A", "B", "C", "D"), "{A,B,C,D}")
        self.assertEqual(self.database.create_postgres_array(1, 2, 3, 4), "{1,2,3,4}")

    def test_add_model_info_and_get_model(self):
        str, param = self.database.add_model_info("testmethod", "testbasis", True)
        self.database.execute(str, param)

        models = self.database.get_models()

        self.assertIn("testmethod/testbasis/True", models)
        self.assertNotIn("wrongmethod/testbasis/True", models)
        self.assertNotIn("testmethod/wrongbasis/True", models)
        self.assertNotIn("testmethod/testbasis/False", models)


        str, param = self.database.add_model_info("testmethod", "testbasis", False)
        self.database.execute(str, param)

        models = self.database.get_models()

        self.assertIn("testmethod/testbasis/True", models)
        self.assertIn("testmethod/testbasis/False", models)

    def test_add_molecule_and_get_molecule(self):

        for i in range(1):

            molecule = self.get_water_monomer()

            str, param = self.database.add_molecule(molecule)

            self.database.execute(str, param)

            self.assertEqual(molecule.get_SHA1(), self.database.get_molecule(molecule.get_SHA1()).get_SHA1())

        for i in range(0):
            molecule = self.get_water_dimer()

            str, param = self.database.add_molecule(molecule)
            self.database.execute(str, param)

            self.assertEqual(molecule, self.database.get_molecule(molecule.get_SHA1()))

    def test_add_calculation_and_get_all_calculations(self):

        calculations = list(self.database.get_all_calculations("testclient", "notused"))
        self.assertEqual(len(calculations), 0)

        molecules = []
        for i in range(100):
            molecules.append(self.get_water_monomer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "tag1")

        calculations = list(self.database.get_all_calculations("testclient", "notused"))
        self.assertEqual(len(calculations), 100)

        for molecule in molecules:
            self.assertIn([molecule, "testmethod", "testbasis", True, False, [0]], calculations)

        calculations = list(self.database.get_all_calculations("testclient", "notused"))
        self.assertEqual(len(calculations), 0)

        molecules = []
        for i in range(100):
            molecules.append(self.get_water_monomer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "tag1")

        calculations = list(self.database.get_all_calculations("testclient", "notused"))
        self.assertEqual(len(calculations), 100)

        for molecule in molecules:
            self.assertIn([molecule, "testmethod", "testbasis", False, False, [0]], calculations)

        calculations = list(self.database.get_all_calculations("testclient", "notused"))
        self.assertEqual(len(calculations), 0)

        molecules = []
        for i in range(100):
            molecules.append(self.get_water_dimer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "tag1")

        calculations = list(self.database.get_all_calculations("testclient", "notused"))
        self.assertEqual(len(calculations), 500)

        for molecule in molecules:
            self.assertIn([molecule, "testmethod", "testbasis", True, False, [0]], calculations)
            self.assertIn([molecule, "testmethod", "testbasis", True, False, [1]], calculations)
            self.assertIn([molecule, "testmethod", "testbasis", True, True, [0]], calculations)
            self.assertIn([molecule, "testmethod", "testbasis", True, True, [1]], calculations)
            self.assertIn([molecule, "testmethod", "testbasis", True, False, [0, 1]], calculations)

        calculations = list(self.database.get_all_calculations("testclient", "notused"))
        self.assertEqual(len(calculations), 0)

        molecules = []
        for i in range(100):
            molecules.append(self.get_water_dimer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", False, "tag1")

        calculations = list(self.database.get_all_calculations("testclient", "notused"))
        self.assertEqual(len(calculations), 300)

        for molecule in molecules:
            self.assertIn([molecule, "testmethod", "testbasis", True, False, [0]], calculations)
            self.assertIn([molecule, "testmethod", "testbasis", True, False, [1]], calculations)
            self.assertIn([molecule, "testmethod", "testbasis", True, False, [0, 1]], calculations)

        calculations = list(self.database.get_all_calculations("testclient", "notused"))
        self.assertEqual(len(calculations), 0)

    def test_build_empty_molecule(self):

        molecule = self.get_water_monomer()

        str, param = self.database.add_molecule(molecule)

        self.database.execute(str, param)

        for atom in molecule.get_atoms():
            atom.set_xyz(0, 0, 0)

        empty_mol = self.database.build_empty_molecule(molecule.get_name())

        self.assertEqual(molecule, empty_mol)

    """
    def test_set_properties_and_get_training_set(self):
        molecules = []
        for i in range(100):
            molecules.append(self.get_water_monomer())

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "tag1")

        calculations = self.database.get_all_calculations("testclient", "notused")

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            calculation_results.append([molecule, method, basis, cp, use_cp, frag_indices, True, random.random(), "some log test"])

        self.database.set_properties(calculation_results)

        training_set = self.database.get_1B_training_set("H2O", "testmethod", "testbasis", True, "tag1")
        self.assertEqual(len(training_set), 100)
        for molecule, method, basis, cp, use_cp, frag_indices, result, energy, log_text in calculation_results:
            self.assertIn([molecule, energy], training_set)


        molecules = []
        energies = []
        for i in range(100):
            molecules.append(self.get_water_dimer())
            energies.append([random.random(), random.random(), random.random()])

        self.database.add_calculations(molecules, "testmethod", "testbasis", True, "tag1")

        calculations = self.database.get_all_calculations("testclient", "notused")

        calculation_results = []

        for molecule, method, basis, cp, use_cp, frag_indices in calculations:
            index = molecules.index(molecule)
            if frag_indices == [0]:
                energy = energies[index][0]
            elif frag_indices == [1]:
                energy = energies[index][1]
            else:
                energy = energies[index][2]
            calculation_results.append([molecule, method, basis, cp, use_cp, frag_indices, True, energy, "some log test"])

        self.database.set_properties(calculation_results)

        training_set = self.database.get_2B_training_set("H2O-H2O", "testmethod", "testbasis", True, "tag1")
        self.assertEqual(len(training_set), 100)

        for molecule, method, basis, cp, use_cp, frag_indices, result, energy, log_text in calculation_results:
            index = molecules.index(molecule)
            self.assertIn([molecule, energy], training_set)
    """


suite = unittest.TestLoader().loadTestsFromTestCase(TestDatabase)