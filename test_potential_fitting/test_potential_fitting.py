import unittest
from . import test_molecule, test_database

suite = unittest.TestSuite([test_molecule.suite, test_database.suite])
