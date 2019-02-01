import unittest
from . import test_model, test_calculator

suite = unittest.TestSuite([test_model.suite, test_calculator.suite])