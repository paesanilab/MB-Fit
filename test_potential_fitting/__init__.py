import unittest
from . import test_potential_fitting

def execute_tests():
    unittest.TextTestRunner(verbosity=2).run(test_potential_fitting.suite)
