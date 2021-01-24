import unittest
from . import test_molecule, test_database, test_calculator, test_utils, test_polynomials, test_fitting, test_case_with_id 

suite = unittest.TestSuite([test_molecule.suite, test_database.suite, test_calculator.suite, test_utils.suite, test_polynomials.suite, test_fitting.suite])
