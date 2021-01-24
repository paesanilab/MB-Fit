import unittest
from .test_mbfit import *

def execute_tests():
    unittest.TextTestRunner(verbosity=2).run(test_mbfit.suite)
