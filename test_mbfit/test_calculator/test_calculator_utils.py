import unittest, os

from test_mbfit.test_case_with_id import TestCaseWithId
from mbfit.calculator import get_calculator, Psi4Calculator, QchemCalculator, fill_energies
from mbfit.exceptions import NoSuchLibraryError

from mbfit.molecule import parse_training_set_file


class TestCalculatorUtils(TestCaseWithId):
    def __init__(self, *args, **kwargs):
        super(TestCalculatorUtils, self).__init__(*args, **kwargs)
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

    def setUpClass():

        TestCalculatorUtils.psi4_settings_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "psi4.ini")
        TestCalculatorUtils.qchem_settings_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "qchem.ini")
        TestCalculatorUtils.bad_code_settings_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "bad_code.ini")

        TestCalculatorUtils.configs_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "configs.xyz")
        TestCalculatorUtils.opt_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2monomer.xyz")
        TestCalculatorUtils.CO2_settings_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "CO2.ini")

        TestCalculatorUtils.training_set_output = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "training_set.xyz")
        TestCalculatorUtils.training_set_reference = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference", "training_set.xyz")

    def test_get_calculator_psi4(self):
        calc = get_calculator(TestCalculatorUtils.psi4_settings_path)

        self.assertIsInstance(calc, Psi4Calculator)

        self.test_passed = True

    def test_get_calculator_qchem(self):

        calc = get_calculator(TestCalculatorUtils.qchem_settings_path)

        self.assertIsInstance(calc, QchemCalculator)

        self.test_passed = True

    def test_get_calculator_bad_code(self):

        with self.assertRaises(NoSuchLibraryError):
            get_calculator(TestCalculatorUtils.bad_code_settings_path)

        self.test_passed = True

    def test_fill_energies(self):
        fill_energies(TestCalculatorUtils.CO2_settings_path,
                      TestCalculatorUtils.configs_path,
                      [TestCalculatorUtils.CO2_settings_path],
                      [TestCalculatorUtils.opt_path],
                      TestCalculatorUtils.training_set_output,
                      "HF",
                      "STO-3G",
                      False)

        out_mols = parse_training_set_file(TestCalculatorUtils.training_set_output)
        ref_mols = parse_training_set_file(TestCalculatorUtils.training_set_reference)

        output = []
        reference = []

        with open(TestCalculatorUtils.training_set_output) as output_file, \
             open(TestCalculatorUtils.training_set_reference) as reference_file:

            out_lines = output_file.readlines()
            ref_lines = reference_file.readlines()

            for i in range(3):
                output.append((round(float(out_lines[i * 5 + 1].split()[0]), 5),
                               round(float(out_lines[i * 5 + 1].split()[1]), 5),
                               out_mols.__next__()))

                reference.append((round(float(ref_lines[i * 5 + 1].split()[0]), 5),
                                  round(float(ref_lines[i * 5 + 1].split()[1]), 5),
                                  ref_mols.__next__()))

        for molecule in output:
            self.assertIn(molecule, reference)
        for molecule in reference:
            self.assertIn(molecule, output)

        self.test_passed = True


suite = unittest.TestLoader().loadTestsFromTestCase(TestCalculatorUtils)
