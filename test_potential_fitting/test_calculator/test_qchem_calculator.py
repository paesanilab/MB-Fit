import unittest
import os

from potential_fitting.calculator import QchemCalculator, Model
from potential_fitting.utils import system
from potential_fitting.exceptions import CommandExecutionError

from .test_calculator import TestCalculator

def hasQchem():

    try:
        system.call("which", "qchem")
    except CommandExecutionError:
        return False
    return True

# this skips this entire test class if qchem is not isntalled.
@unittest.skipUnless(hasQchem(), "Qchem is not installed and cannot be tested.")
class TestQchemCalculator(TestCalculator):

    # set up before each test case
    def setUp(self):
        self.calculator1 = QchemCalculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), True)
        self.calculator2 = QchemCalculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "SCN-monomer.ini"), True)
        self.calculator3 = QchemCalculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), True)
        self.calculator4 = QchemCalculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "SCN-monomer.ini"), True)

    def test_calculate_energy(self):

        
        self.calculator1.calculate_energy(TestCalculator.CO2, TestCalculator.model1, [0])

        self.calculator2.calculate_energy(TestCalculator.SCN, TestCalculator.model1, [0])

        self.calculator3.calculate_energy(TestCalculator.CO2, TestCalculator.model2, [0])

        self.calculator4.calculate_energy(TestCalculator.SCN, TestCalculator.model2, [0])

    def test_optimize_geometry(self):
        
        self.calculator1.optimize_geometry(TestCalculator.CO2, TestCalculator.model1)

        self.calculator2.optimize_geometry(TestCalculator.SCN, TestCalculator.model1)

        self.calculator3.optimize_geometry(TestCalculator.CO2, TestCalculator.model2)

        self.calculator4.optimize_geometry(TestCalculator.SCN, TestCalculator.model2)

    def test_find_frequencies(self):
        
        self.calculator1.find_frequencies(TestCalculator.CO2, TestCalculator.model1)

        self.calculator2.find_frequencies(TestCalculator.SCN, TestCalculator.model1)

        self.calculator3.find_frequencies(TestCalculator.CO2, TestCalculator.model2)

        self.calculator4.find_frequencies(TestCalculator.SCN, TestCalculator.model2)

suite = unittest.TestLoader().loadTestsFromTestCase(TestQchemCalculator)