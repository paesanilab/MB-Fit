import unittest
import os

class TestCaseWithId(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestCaseWithId, self).__init__(*args, **kwargs)
        self.test_passed = False
        self.test_name = self.id()
        self.test_folder = ""

    # clean up after each test case
    def tearDown(self):
        local_output = os.path.join(self.test_folder,"output")
        mbfithome = os.environ.get('MBFIT_HOME')
        if os.path.isdir(local_output):
            if self.test_passed:
                os.system("mkdir -p " + os.path.join(mbfithome, "passed_tests_outputs"))
                os.system("mv " + local_output + " " + os.path.join(mbfithome, "passed_tests_outputs", self.test_name))
            else:
                os.system("mkdir -p " + os.path.join(mbfithome, "failed_tests_outputs"))
                os.system("mv " + local_output + " " + os.path.join(mbfithome, "failed_tests_outputs", self.test_name))

