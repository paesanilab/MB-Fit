import unittest
import os


from potential_fitting.molecule import Molecule
from potential_fitting.calculator import Calculator, Model

class TestCalculator(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestCalculator, self).__init__(*args, **kwargs)
        self.test_passed = False
        self.test_name = self.id()

    # set up before the first test case
    def setUpClass():

        TestCalculator.CO2 = Molecule.read_xyz_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.xyz"), [3], ["CO2"], [0], [1], ["A1B2"], ["C(O)O"])
        TestCalculator.CN = Molecule.read_xyz_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.xyz"), [2], ["CN"], [-1], [1], ["A1B1"], ["CN"])

        TestCalculator.model1 = Model("HF", "STO-3G", True)
        TestCalculator.model2 = Model("wb97", "cc-pvdz", False)

    # clean up after the last test case
    def tearDownClass():
        pass

    # set up before each test case
    def setUp(self):
        self.calculator1 = Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), False)
        self.calculator2 = Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.ini"), False)
        self.calculator3 = Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), False)
        self.calculator4 = Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.ini"), False)

    # clean up after each test case
    def tearDown(self):
        local_output = os.path.join(os.path.dirname(os.path.abspath(__file__)),"output")
        mbfithome = os.environ.get('MBFIT_HOME')
        if os.path.isdir(local_output):
            if self.test_passed:
                os.system("mkdir -p " + os.path.join(mbfithome, "passed_tests_outputs"))
                os.system("mv " + local_output + " " + os.path.join(mbfithome, "passed_tests_outputs", self.test_name))
            else:
                os.system("mkdir -p " + os.path.join(mbfithome, "failed_tests_outputs"))
                os.system("mv " + local_output + " " + os.path.join(mbfithome, "failed_tests_outputs", self.test_name))

    def test_set_logging(self):
        self.test_passed = True

    def test_is_installed(self):
        with self.assertRaises(NotImplementedError):
            self.calculator1.is_installed()

        with self.assertRaises(NotImplementedError):
            self.calculator2.is_installed()

        with self.assertRaises(NotImplementedError):
            self.calculator2.is_installed()

        with self.assertRaises(NotImplementedError):
            self.calculator2.is_installed()

        self.test_passed = True

    def test_calculate_energy(self):

        with self.assertRaises(NotImplementedError):
            self.calculator1.calculate_energy(TestCalculator.CO2, TestCalculator.model1, [0])

        with self.assertRaises(NotImplementedError):
            self.calculator2.calculate_energy(TestCalculator.CN, TestCalculator.model1, [0])

        with self.assertRaises(NotImplementedError):
            self.calculator3.calculate_energy(TestCalculator.CO2, TestCalculator.model2, [0])

        with self.assertRaises(NotImplementedError):
            self.calculator4.calculate_energy(TestCalculator.CN, TestCalculator.model2, [0])

        self.test_passed = True

    def test_optimize_geometry(self):
        
        with self.assertRaises(NotImplementedError):
            self.calculator1.optimize_geometry(TestCalculator.CO2, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator2.optimize_geometry(TestCalculator.CN, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator3.optimize_geometry(TestCalculator.CO2, TestCalculator.model2)

        with self.assertRaises(NotImplementedError):
            self.calculator4.optimize_geometry(TestCalculator.CN, TestCalculator.model2)

        self.test_passed = True

    def test_calculate_frequencies(self):
        
        with self.assertRaises(NotImplementedError):
            self.calculator1.calculate_frequencies(TestCalculator.CO2, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator2.calculate_frequencies(TestCalculator.CN, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator3.calculate_frequencies(TestCalculator.CO2, TestCalculator.model2)

        with self.assertRaises(NotImplementedError):
            self.calculator4.calculate_frequencies(TestCalculator.CN, TestCalculator.model2)

        self.test_passed = True

    def test_is_valid_model(self):

        with self.assertRaises(NotImplementedError):
            self.calculator1.is_valid_model(Model("HF", "STO-3G", False))

        with self.assertRaises(NotImplementedError):
            self.calculator2.is_valid_model(Model("wb97m-v", "STO-3G", True))

        with self.assertRaises(NotImplementedError):
            self.calculator3.is_valid_model(Model("", "", False))

        with self.assertRaises(NotImplementedError):
            self.calculator4.is_valid_model(Model("abcd", "efgh", False))

        self.test_passed = True

    def test_check_neg_freqs(self):
        
        self.assertEquals(self.calculator1.check_neg_freqs([]), 0)
        self.assertEquals(self.calculator1.check_neg_freqs([0, 1, 2, 3, 4]), 0)
        self.assertEquals(self.calculator1.check_neg_freqs([0, 0, -1, 1]), 1)
        self.assertEquals(self.calculator1.check_neg_freqs([-1, -2, -3, -0.4, -234234]), 5)
        self.assertEquals(self.calculator1.check_neg_freqs([-2, -2, 2, 2, 4, 1234567890, 0.1234567890, -0.0000001]), 3)
        self.test_passed = True

suite = unittest.TestLoader().loadTestsFromTestCase(TestCalculator)
