import unittest
import testMolecule
import sys

# declare test suite of all tests
suite = unittest.TestSuite([testMolecule.suite])

# execute the test suite
result = unittest.TextTestRunner(verbosity=2).run(suite)

# check if any failures or errors occured
if len(result.failures) != 0 or len(result.errors) != 0:
    test_log = open("tests.log", "w+") # w+ opens file to write and creates if it didn't already exist
    for failure in result.failures:
        test_log.write(failure)
    for error in result.errors:
        test_log.write(error)
    test_log.close();

    sys.exit(1) # exit with failure exit code
else:
    test_log = open("tests.log", "w+") # w+ opens file to write and creates if it didn't already exist
    test_log.write("ALL TESTS PASSED! WOO! :DDD \n")
    test_log.close()

    sys.exit(0) # exit with success exit code


