import unittest
from . import test_atom, test_fragment, test_molecule

suite = unittest.TestSuite([test_atom.suite, test_fragment.suite, test_molecule.suite])
