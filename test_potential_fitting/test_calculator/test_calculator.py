import unittest
import os


from potential_fitting.molecule import Molecule
from potential_fitting.calculator import Calculator, Model

class TestCalculator(unittest.TestCase):

    # set up before the first test case
    def setUpClass():

        TestCalculator.CO2 = Molecule().read_xyz_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.xyz"), [3], ["CO2"], [0], [1], ["A1B2"])
        TestCalculator.CN = Molecule().read_xyz_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.xyz"), [2], ["CN"], [-1], [1], ["A1B1"])

        TestCalculator.model1 = Model("HF", "STO-3G", True)
        TestCalculator.model2 = Model("wb97", "cc-pvdz", False)
        pass

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
        pass

    def test_set_logging(self):
        pass

    def test_calculate_energy(self):

        with self.assertRaises(NotImplementedError):
            self.calculator1.calculate_energy(TestCalculator.CO2, TestCalculator.model1, [0])

        with self.assertRaises(NotImplementedError):
            self.calculator2.calculate_energy(TestCalculator.CN, TestCalculator.model1, [0])

        with self.assertRaises(NotImplementedError):
            self.calculator3.calculate_energy(TestCalculator.CO2, TestCalculator.model2, [0])

        with self.assertRaises(NotImplementedError):
            self.calculator4.calculate_energy(TestCalculator.CN, TestCalculator.model2, [0])

    def test_optimize_geometry(self):
        
        with self.assertRaises(NotImplementedError):
            self.calculator1.optimize_geometry(TestCalculator.CO2, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator2.optimize_geometry(TestCalculator.CN, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator3.optimize_geometry(TestCalculator.CO2, TestCalculator.model2)

        with self.assertRaises(NotImplementedError):
            self.calculator4.optimize_geometry(TestCalculator.CN, TestCalculator.model2)

    def test_find_frequencies(self):
        
        with self.assertRaises(NotImplementedError):
            self.calculator1.find_frequencies(TestCalculator.CO2, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator2.find_frequencies(TestCalculator.CN, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator3.find_frequencies(TestCalculator.CO2, TestCalculator.model2)

        with self.assertRaises(NotImplementedError):
            self.calculator4.find_frequencies(TestCalculator.CN, TestCalculator.model2)

    def test_check_neg_freqs(self):
        
        self.assertEquals(self.calculator1.check_neg_freqs([]), 0)
        self.assertEquals(self.calculator1.check_neg_freqs([0, 1, 2, 3, 4]), 0)
        self.assertEquals(self.calculator1.check_neg_freqs([0, 0, -1, 1]), 1)
        self.assertEquals(self.calculator1.check_neg_freqs([-1, -2, -3, -0.4, -234234]), 5)
        self.assertEquals(self.calculator1.check_neg_freqs([-2, -2, 2, 2, 4, 1234567890, 0.1234567890, -0.0000001]), 3)

suite = unittest.TestLoader().loadTestsFromTestCase(TestCalculator)