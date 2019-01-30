import unittest

from potential_fitting.calculator import Calcualtor, Model

class TestCalculator(unittest.TestCase):

    # set up before the first test case
    def setUpClass():

        TestCalculator.H2O = None # construct the molecule
        TestCalculator.NO3 = None # Constructu the molecule

        TestCalculator.model1 = Model("HF", "STO-3G", True)
        TestCalculator.model2 = Model("wb97m-v", "cc-pvdz", False)
        pass

    # clean up after the last test case
    def tearDownClass():
        pass

    # set up before each test case
    def setUp(self):
        self.calculator1 = Calculator("settings1.ini", False)
        self.calculator2 = Calculator("settings2.ini". False)
        self.calculator3 = Calculator("settings2.ini". False)
        self.calculator4 = Calculator("settings2.ini". False)

    # clean up after each test case
    def tearDown(self):
        pass

    def test_set_logging(self):
        pass

    def test_calculate_energy(self):

        with self.assertRaises(NotImplementedError):
            self.calculator1.calculate_energy(TestCalculator.H2O, TestCalculator.model1, [0])

        with self.assertRaises(NotImplementedError):
            self.calculator2.calculate_energy(TestCalculator.NO3, TestCalculator.model1, [])

        with self.assertRaises(NotImplementedError):
            self.calculator3.calculate_energy(TestCalculator.H2O, TestCalculator.model2, [])

        with self.assertRaises(NotImplementedError):
            self.calculator4.calculate_energy(TestCalculator.NO3, TestCalculator.model2, [0])

    def test_optimize_geometry(self):
        
        with self.assertRaises(NotImplementedError):
            self.calculator1.optimize_geometry(TestCalculator.H2O, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator2.optimize_geometry(TestCalculator.NO3, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator3.optimize_geometry(TestCalculator.H2O, TestCalculator.model2)

        with self.assertRaises(NotImplementedError):
            self.calculator4.optimize_geometry(TestCalculator.NO3, TestCalculator.model2)

    def test_find_frequencies(self):
        
        with self.assertRaises(NotImplementedError):
            self.calculator1.find_frequencies(TestCalculator.H2O, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator2.find_frequencies(TestCalculator.NO3, TestCalculator.model1)

        with self.assertRaises(NotImplementedError):
            self.calculator3.find_frequencies(TestCalculator.H2O, TestCalculator.model2)

        with self.assertRaises(NotImplementedError):
            self.calculator4.find_frequencies(TestCalculator.NO3, TestCalculator.model2)

    def test_check_neg_freqs(self):
        
        self.assertEquals(self.calculator1.check_neg_freqs([]), 0)
        self.assertEquals(self.calculator1.check_neg_freqs([0, 1, 2, 3, 4]), 0)
        self.assertEquals(self.calculator1.check_neg_freqs([0, 0, -1, 1]), -1)
        self.assertEquals(self.calculator1.check_neg_freqs([-1, -2, -3, -0.4, -234234]), 5)
        self.assertEquals(self.calculator1.check_neg_freqs([-2, -2, 2, 2, 4, 1234567890, 0.1234567890, -0.0000001]), 3)

suite = unittest.TestLoader().loadTestsFromTestCase(TestCalculator)