import unittest

from potential_fitting.calculator import Psi4Calcualtor, Model

from .test_calculator import TestCalculator

class TestPsi4Calculator(TestCalculator):

    # set up before each test case
    def setUp(self):
        self.calculator1 = Psi4Calcualtor("settings1.ini", False)
        self.calculator2 = Psi4Calcualtor("settings2.ini". False)
        self.calculator3 = Psi4Calcualtor("settings2.ini". False)
        self.calculator4 = Psi4Calcualtor("settings2.ini". False)

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

suite = unittest.TestLoader().loadTestsFromTestCase(TestCalculator)