import unittest, os

from test_potential_fitting.test_case_with_id import TestCaseWithId
from potential_fitting.utils import math

class TestMath(TestCaseWithId):

    def test_test_difference_under_threshold(self):
        self.assertTrue(math.test_difference_under_threshold(1, 2, 1.1))
        self.assertFalse(math.test_difference_under_threshold(1, 2, 1))


        self.assertTrue(math.test_difference_under_threshold(1, 2, 2))
        self.assertFalse(math.test_difference_under_threshold(1, 4, 2))

        self.assertTrue(math.test_difference_under_threshold(1.00001, 1.00002, 0.00002))
        self.assertFalse(math.test_difference_under_threshold(1.00001, 1.00002, 0.000009)) # uses 0.000009 instead of 0.00001 due to float innacuracy.

        self.assertTrue(math.test_difference_under_threshold(1, 1, 0.0002))
        self.assertFalse(math.test_difference_under_threshold(1, 1, 0))

        self.assertTrue(math.test_difference_under_threshold(1, 1.1, 5))
        self.assertFalse(math.test_difference_under_threshold(1, 5.1, 4))

        self.test_passed = True

suite = unittest.TestLoader().loadTestsFromTestCase(TestMath)
