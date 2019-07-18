import unittest
from . import test_generate_input_poly, test_generate_poly

suite = unittest.TestSuite([test_generate_input_poly.suite, test_generate_poly.suite])
