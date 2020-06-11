import unittest, os

from potential_fitting.utils import system, files
from potential_fitting.exceptions import CommandExecutionError, CommandNotFoundError

class TestSystem(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestSystem, self).__init__(*args, **kwargs)
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

    def setUpClass():
        TestSystem.bad_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "bad_file.txt")
        TestSystem.in_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "in_file.txt")
        TestSystem.out_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "out_file.txt")

    def test_call(self):

        files.init_file(TestSystem.in_file)
        files.init_file(TestSystem.out_file)

        with self.assertRaises(CommandExecutionError):
            system.call("cat", TestSystem.bad_file)

        with self.assertRaises(CommandNotFoundError):
            system.call("bad_command")

        with open(TestSystem.out_file, "w") as out_file:
            system.call("echo", "echo echo  echo   echo   echo    echo", out_file=out_file)

        with open(TestSystem.out_file, "r") as out_file:
            self.assertEqual(out_file.read(), "echo echo  echo   echo   echo    echo\n")

        with open(TestSystem.in_file, "w") as in_file:
            in_file.write("cat!\n")

        with open(TestSystem.in_file, "r") as in_file, open(TestSystem.out_file, "w") as out_file:
            system.call("cat", in_file=in_file, out_file=out_file)

        with open(TestSystem.out_file, "r") as out_file:
            self.assertEqual(out_file.read(), "cat!\n")

        self.test_passed = True

suite = unittest.TestLoader().loadTestsFromTestCase(TestSystem)
