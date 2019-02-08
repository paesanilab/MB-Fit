import unittest
import os

from potential_fitting.calculator import Psi4Calculator, Model

from .test_calculator import TestCalculator

class TestPsi4Calculator(TestCalculator):

    # set up before each test case
    def setUp(self):
        self.calculator1 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), True)
        self.calculator2 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.ini"), True)
        self.calculator3 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), True)
        self.calculator4 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.ini"), True)

    def test_calculate_energy(self):
        
        self.calculator1.calculate_energy(TestCalculator.CO2, TestCalculator.model1, [0])

        self.calculator2.calculate_energy(TestCalculator.CN, TestCalculator.model1, [0])

        self.calculator3.calculate_energy(TestCalculator.CO2, TestCalculator.model2, [0])

        self.calculator4.calculate_energy(TestCalculator.CN, TestCalculator.model2, [0])

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