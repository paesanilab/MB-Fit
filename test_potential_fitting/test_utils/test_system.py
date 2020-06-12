import unittest, os

from test_potential_fitting.test_case_with_id import TestCaseWithId
from potential_fitting.utils import system, files
from potential_fitting.exceptions import CommandExecutionError, CommandNotFoundError

class TestSystem(TestCaseWithId):

    def setUpClass():
        TestSystem.bad_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "bad_file.txt")
        TestSystem.in_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "in_file.txt")
        TestSystem.out_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "out_file.txt")

    def test_call(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

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
