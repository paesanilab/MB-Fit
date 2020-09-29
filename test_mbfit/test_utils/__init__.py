import unittest
from . import test_constants, test_math, test_files, test_system, test_settings_reader, test_quaternion

suite = unittest.TestSuite([test_constants.suite, test_math.suite, test_files.suite, test_system.suite, test_settings_reader.suite, test_quaternion.suite])
