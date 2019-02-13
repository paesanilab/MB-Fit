import unittest
import os

from potential_fitting.calculator import Psi4Calculator, Model
from potential_fitting.utils import math, files

from .test_calculator import TestCalculator

def hasPsi4():

    try:
        import psi4
    except ImportError:
        return False
    return True

@unittest.skipUnless(hasPsi4(), "Psi4 is not installed and cannot be tested.")
class TestPsi4Calculator(TestCalculator):

    # set up before each test case
    def setUp(self):
        self.calculator1 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), True)
        self.calculator2 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.ini"), True)
        self.calculator3 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), True)
        self.calculator4 = Psi4Calculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.ini"), True)

    def test_is_installed(self):
        try:
            import psi4
            self.assertTrue(self.calculator1.is_installed())
            self.assertTrue(self.calculator2.is_installed())
            self.assertTrue(self.calculator3.is_installed())
            self.assertTrue(self.calculator4.is_installed())
        except ImportError:
            self.assertFalse(self.calculator1.is_installed())
            self.assertFalse(self.calculator2.is_installed())
            self.assertFalse(self.calculator3.is_installed())
            self.assertFalse(self.calculator4.is_installed())

    def test_calculate_energy(self):
        energy, log_path = self.calculator1.calculate_energy(TestCalculator.CO2, TestCalculator.model1, [0])

        ref_energy, ref_log_path = -184.825948265526, "869cbba7_2019-02-12_17-02-55.820866.out"

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Energy calculation failed. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

        energy, log_path = self.calculator2.calculate_energy(TestCalculator.CN, TestCalculator.model1, [0])

        ref_energy, ref_log_path = -90.8324300259070, "36c6e9c4_2019-02-12_17-27-48.298430.out"

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Energy calculation failed. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

        energy, log_path = self.calculator3.calculate_energy(TestCalculator.CO2, TestCalculator.model2, [0])

        ref_energy, ref_log_path = -188.398728493588, "869cbba7_2019-02-12_16-57-06.189857.out"

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Energy calculation failed. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

        energy, log_path = self.calculator4.calculate_energy(TestCalculator.CN, TestCalculator.model2, [0])

        ref_energy, ref_log_path = -92.7149022432185, "36c6e9c4_2019-02-12_17-25-30.184407.out"

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Energy calculation failed. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

    def test_optimize_geometry(self):

        # compare energy and geomtery in standard orientation
        
        self.calculator1.optimize_geometry(TestCalculator.CO2, TestCalculator.model1)

        self.calculator2.optimize_geometry(TestCalculator.CN, TestCalculator.model1)

        self.calculator3.optimize_geometry(TestCalculator.CO2, TestCalculator.model2)

        self.calculator4.optimize_geometry(TestCalculator.CN, TestCalculator.model2)

    def test_find_frequencies(self):

        # compare freqs, normal modes, and reduced masses.
        
        self.calculator1.find_frequencies(TestCalculator.CO2, TestCalculator.model1)

        self.calculator2.find_frequencies(TestCalculator.CN, TestCalculator.model1)

        self.calculator3.find_frequencies(TestCalculator.CO2, TestCalculator.model2)

        self.calculator4.find_frequencies(TestCalculator.CN, TestCalculator.model2)

suite = unittest.TestLoader().loadTestsFromTestCase(TestPsi4Calculator)