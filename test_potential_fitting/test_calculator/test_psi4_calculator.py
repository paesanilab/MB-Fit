import unittest
import os

from potential_fitting.calculator import Psi4Calculator, Model
from potential_fitting.utils import math, files

from .test_calculator import TestCalculator

class TestPsi4Calculator(TestCalculator):

    # set up before each test case
    def setUp(self):
        self.calculator1 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), True)
        self.calculator2 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.ini"), True)
        self.calculator3 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), True)
        self.calculator4 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.ini"), True)

    def test_calculate_energy(self):
        
        energy, log_path = self.calculator1.calculate_energy(TestCalculator.CO2, TestCalculator.model1, [0])

        ref_energy = -184.825948265526

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1/10000))

        energy, log_path = self.calculator2.calculate_energy(TestCalculator.CN, TestCalculator.model1, [0])

        ref_energy = -90.8324300259072

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1/10000))

        energy, log_path = self.calculator3.calculate_energy(TestCalculator.CO2, TestCalculator.model2, [0])

        ref_energy = -188.398728493588

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1/10000))

        energy, log_path = self.calculator4.calculate_energy(TestCalculator.CN, TestCalculator.model2, [0])

        ref_energy = -92.7149022432186

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1/10000))

    def test_optimize_geometry(self):
        
        self.calculator1.optimize_geometry(TestCalculator.CO2, TestCalculator.model1)

        self.calculator2.optimize_geometry(TestCalculator.CN, TestCalculator.model1)

        self.calculator3.optimize_geometry(TestCalculator.CO2, TestCalculator.model2)

        self.calculator4.optimize_geometry(TestCalculator.CN, TestCalculator.model2)

    def test_find_frequencies(self):
        
        self.calculator1.find_frequencies(TestCalculator.CO2, TestCalculator.model1)

        self.calculator2.find_frequencies(TestCalculator.CN, TestCalculator.model1)

        self.calculator3.find_frequencies(TestCalculator.CO2, TestCalculator.model2)

        self.calculator4.find_frequencies(TestCalculator.CN, TestCalculator.model2)

suite = unittest.TestLoader().loadTestsFromTestCase(TestPsi4Calculator)