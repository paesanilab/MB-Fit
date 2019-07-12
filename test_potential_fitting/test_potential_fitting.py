import unittest
from . import test_molecule, test_database, test_calculator, test_utils

suite = unittest.TestSuite([test_molecule.suite, test_database.suite, test_calculator.suite, test_utils.suite])
