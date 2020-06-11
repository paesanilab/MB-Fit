import unittest, os

from potential_fitting.utils import math

class TestMath(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestMath, self).__init__(*args, **kwargs)
        self.test_passed = False
        self.test_name = self.id()

    # clean up after each test case
    def tearDown(self):
        local_output = os.path.join(os.path.dirname(os.path.abspath(__file__)),"output")
        mbfithome = os.environ.get('MBFIT_HOME')
        if os.path.isdir(local_output):
            if self.test_passed:
                os.system("mkdir -p " + os.path.join(mbfithome, "passed_tests_outputs"))
                os.system("mv " + local_output + " " + os.path.join(mbfithome, "passed_tests_outputs", self.test_name))
            else:
                os.system("mkdir -p " + os.path.join(mbfithome, "failed_tests_outputs"))
                os.system("mv " + local_output + " " + os.path.join(mbfithome, "failed_tests_outputs", self.test_name))

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
