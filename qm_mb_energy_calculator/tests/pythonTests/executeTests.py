import unittest
import testMolecule
import sys

# declare test suite of all tests
suite = unittest.TestSuite([testMolecule.suite])

# execute the test suite
result = unittest.TextTestRunner(verbosity=2).run(suite)

if len(result.failures) != 0 or len(result.errors) != 0:
    sys.exit(1)
else:
    sys.exit(0)


