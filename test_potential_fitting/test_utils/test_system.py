import unittest

from potential_fitting.utils import system
from potential_fitting.exceptions import CommandExecutionError, CommandNotFoundError

class TestSystem(unittest.TestCase):

    def test_call(self):

        with self.assertRaises(CommandExecutionError):
            system.call("cat", "bad_file.txt")

        with self.assertRaises(CommandNotFoundError):
            system.call("bad_command")

        with open("out_file.txt", "w") as out_file:
            system.call("echo", "echo echo  echo   echo   echo    echo", out_file=out_file)

        with open("out_file.txt", "r") as out_file:
            self.assertEqual(out_file.read(), "echo echo  echo   echo   echo    echo\n")

        with open("in_file.txt", "w") as in_file:
            in_file.write("cat!\n")

        with open("in_file.txt", "r") as in_file, open("out_file.txt", "w") as out_file:
            system.call("cat", in_file=in_file, out_file=out_file)

        with open("out_file.txt", "r") as out_file:
            self.assertEqual(out_file.read(), "cat!\n")


suite = unittest.TestLoader().loadTestsFromTestCase(TestSystem)