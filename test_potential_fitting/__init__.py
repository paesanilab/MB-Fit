import unittest
from .test_potential_fitting import *

def execute_tests():
    unittest.TextTestRunner(verbosity=2).run(test_potential_fitting.suite)
