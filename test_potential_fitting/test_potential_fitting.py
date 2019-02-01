import unittest
from . import test_molecule, test_calculator

suite = unittest.TestSuite([test_molecule.suite, test_calculator.suite])
