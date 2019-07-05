import unittest
import os

from potential_fitting.calculator import QchemCalculator, Model
from potential_fitting.utils import system, math
from potential_fitting.exceptions import CommandExecutionError
from potential_fitting.molecule import Molecule

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
        self.calculator2 = QchemCalculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.ini"), True)
        self.calculator3 = QchemCalculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.ini"), True)
        self.calculator4 = QchemCalculator(os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CN-monomer.ini"), True)

    def test_is_installed(self):
        try:
            system.call("which", "qchem")
            self.assertTrue(self.calculator1.is_installed())
            self.assertTrue(self.calculator2.is_installed())
            self.assertTrue(self.calculator3.is_installed())
            self.assertTrue(self.calculator4.is_installed())
        except CommandExecutionError:
            self.assertFalse(self.calculator1.is_installed())
            self.assertFalse(self.calculator2.is_installed())
            self.assertFalse(self.calculator3.is_installed())
            self.assertFalse(self.calculator4.is_installed())

    def test_calculate_energy(self):

        # include log files from reference calculations
        
        energy, log_path = self.calculator1.calculate_energy(TestCalculator.CO2, TestCalculator.model1, [0])

        ref_energy, ref_log_path = -184.8257490397, "reference/869cbba7_2019-02-12_16-50-54.163398.out"

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Energy calculation failed. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

        energy, log_path = self.calculator2.calculate_energy(TestCalculator.CN, TestCalculator.model1, [0])

        ref_energy, ref_log_path = -90.8322075633, "reference/36c6e9c4_2019-02-12_17-14-46.767262.out"

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Energy calculation failed. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

        energy, log_path = self.calculator3.calculate_energy(TestCalculator.CO2, TestCalculator.model2, [0])

        ref_energy, ref_log_path = -188.3943871832, "reference/869cbba7_2019-02-12_16-54-55.906360.out"

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Energy calculation failed. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

        energy, log_path = self.calculator4.calculate_energy(TestCalculator.CN, TestCalculator.model2, [0])

        ref_energy, ref_log_path = -92.7130711210, "reference/36c6e9c4_2019-02-12_17-18-21.567344.out"

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Energy calculation failed. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

    def test_optimize_geometry(self):
        
        geometry, energy, log_path = self.calculator1.optimize_geometry(TestCalculator.CO2, TestCalculator.model1)

        ref_geometry, ref_energy, ref_log_path = Molecule.read_xyz(
                "3\n" +
                "comment line\n" + 
                "C    0.00000000000000e+00   0.00000000000000e+00   0.00000000000000e+00\n" + 
                "O    0.00000000000000e+00   0.00000000000000e+00   1.18790715320000e+00\n" + 
                "O    0.00000000000000e+00   0.00000000000000e+00  -1.18790715320000e+00",
                [3], ["CO2"], [0], [1], ["A1B2"], ["C(O)O"]
                ), -185.068390561172, "reference/87d98bce_2019-02-12_17-36-57.640340.out"

        self.assertTrue(geometry == ref_geometry, "Geometry optimization failed. Compare log files reference: {} and calculated: {}.".format(ref_log_path, log_path))

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Geometry optimization returned incorrect energy. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

        geometry, energy, log_path = self.calculator2.optimize_geometry(TestCalculator.CN, TestCalculator.model1)

        ref_geometry, ref_energy, ref_log_path = Molecule.read_xyz(
                "2\n" + 
                "comment line\n" + 
                "C    0.00000000000000e+00   0.00000000000000e+00  -6.19578775600000e-01\n" + 
                "N    0.00000000000000e+00   0.00000000000000e+00   5.42655698700000e-01\n",
                [2], ["CN"], [-1], [1], ["A1B1"], ["CN"]
                ), -90.9376700608474, "reference/59bb6335_2019-02-26_16-43-16.263184.out"

        self.assertTrue(geometry == ref_geometry, "Geometry optimization failed. Compare log files reference: {} and calculated: {}.".format(ref_log_path, log_path))

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Geometry optimization returned incorrect energy. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

        geometry, energy, log_path = self.calculator3.optimize_geometry(TestCalculator.CO2, TestCalculator.model2)

        ref_geometry, ref_energy, ref_log_path = Molecule.read_xyz(
                "3\n" +
                "comment line\n" + 
                "C    0.00000000000000e+00   0.00000000000000e+00   0.00000000000000e+00\n" + 
                "O    0.00000000000000e+00   0.00000000000000e+00   1.16576119960000e+00\n" + 
                "O    0.00000000000000e+00   0.00000000000000e+00  -1.16576119960000e+00",
                [3], ["CO2"], [0], [1], ["A1B2"], ["C(O)O"]
                ), -188.565745055677, "reference/87d98bce_2019-02-12_17-54-18.734011.out"

        self.assertTrue(geometry == ref_geometry, "Geometry optimization failed. Compare log files reference: {} and calculated: {}.".format(ref_log_path, log_path))

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Geometry optimization returned incorrect energy. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

        geometry, energy, log_path = self.calculator4.optimize_geometry(TestCalculator.CN, TestCalculator.model2)

        ref_geometry, ref_energy, ref_log_path = Molecule.read_xyz(
                "2\n" +
                "comment line\n" + 
                "C    0.00000000000000e+00   0.00000000000000e+00  -6.30033057000000e-01\n" + 
                "N    0.00000000000000e+00   0.00000000000000e+00   5.53109980100000e-01\n" + 
                [2], ["CN"], [-1], [1], ["A1B1"], ["CN"]
                ), -92.8181450462365, "reference/59bb6335_2019-02-26_16-47-33.759389.out"

        self.assertTrue(geometry == ref_geometry, "Geometry optimization failed. Compare log files reference: {} and calculated: {}.".format(ref_log_path, log_path))

        self.assertTrue(math.test_difference_under_threshold(energy, ref_energy, 1e-5), "Geometry optimization returned incorrect energy. Reference: {}, calculated: {}. Compare log files reference: {} and calculated: {}.".format(ref_energy, energy, ref_log_path, log_path))

    def test_find_frequencies(self):
        
        self.calculator1.find_frequencies(TestCalculator.CO2, TestCalculator.model1)

        self.calculator2.find_frequencies(TestCalculator.CN, TestCalculator.model1)

        self.calculator3.find_frequencies(TestCalculator.CO2, TestCalculator.model2)

        self.calculator4.find_frequencies(TestCalculator.CN, TestCalculator.model2)

suite = unittest.TestLoader().loadTestsFromTestCase(TestQchemCalculator)