import unittest
import testMolecule
import sys
sys.path.insert(0, "./monomers")
import testMonomers

# declare test suite of all tests
suite = unittest.TestSuite([testMolecule.suite, testMonomers.suite])

# execute the test suite
result = unittest.TextTestRunner(verbosity=2).run(suite)

# check if any failures or errors occured
if len(result.failures) != 0 or len(result.errors) != 0:
    test_log = open("tests.log", "w+") # w+ opens file to write and creates if it didn't already exist
    if len(result.failures)  != 0:
        test_log.write("TEST FAILURES: \n\n")
        for failure in result.failures:
            test_log.write("TEST {} FAILED:\n".format(failure[0]))
            test_log.write("Detailed Failure Information: \n")
            test_log.write(failure[1] + "\n")
        test_log.write("\n")
    if len(result.errors)  != 0:
        test_log.write("TEST ERRORS: \n\n")
        for error in result.errors:
            test_log.write("TEST {} ERRORED:\n".format(error[0]))
            test_log.write("Detailed Error Information: \n")
            test_log.write(error[1] + "\n")
    test_log.close();

    sys.exit(1) # exit with failure exit code
else:
    test_log = open("tests.log", "w+") # w+ opens file to write and creates if it didn't already exist
    test_log.write("ALL TESTS PASSED! WOO! :DDD \n")
    test_log.close()

    sys.exit(0) # exit with success exit code


