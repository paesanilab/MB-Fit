import unittest
import os

from potential_fitting.polynomials import generate_poly

class TestGeneratePoly(unittest.TestCase):

    def setUpClass():
        TestGeneratePoly.settings = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources", "all.ini")

    def test_A1B2(self):
        pass

suite = unittest.TestLoader().loadTestsFromTestCase(TestGeneratePoly)