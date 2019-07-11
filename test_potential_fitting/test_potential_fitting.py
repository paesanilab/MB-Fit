import unittest
from . import test_molecule, test_database, test_calculator

suite = unittest.TestSuite([test_molecule.suite, test_database.suite, test_calculator.suite])
suite = unittest.TestSuite([test_database.suite])
