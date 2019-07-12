import unittest
from . import test_model, test_calculator, test_psi4_calculator, test_qchem_calculator

suite = unittest.TestSuite([test_model.suite, test_calculator.suite, test_psi4_calculator.suite, test_qchem_calculator.suite])