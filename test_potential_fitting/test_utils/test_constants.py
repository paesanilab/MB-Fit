import unittest

from potential_fitting.utils import constants
from potential_fitting.exceptions import InvalidValueError

class TestConstants(unittest.TestCase):

    def test_symbol_to_number(self):
        self.assertEqual(constants.symbol_to_number("H"), 1)
        self.assertEqual(constants.symbol_to_number("He"), 2)
        self.assertEqual(constants.symbol_to_number("O"), 8)
        self.assertEqual(constants.symbol_to_number("P"), 15)
        self.assertEqual(constants.symbol_to_number("H"), 1)
        self.assertEqual(constants.symbol_to_number("Ne"), 10)

        self.assertEqual(constants.symbol_to_number("NE"), 10)
        self.assertEqual(constants.symbol_to_number("ne"), 10)

        with self.assertRaises(InvalidValueError):
            constants.symbol_to_number("Yo")

        with self.assertRaises(InvalidValueError):
            constants.symbol_to_number("Xx")


    def test_number_to_symbol(self):
        self.assertEqual(constants.number_to_symbol(1), "H")
        self.assertEqual(constants.number_to_symbol(2), "He")
        self.assertEqual(constants.number_to_symbol(8), "O")
        self.assertEqual(constants.number_to_symbol(15), "P")
        self.assertEqual(constants.number_to_symbol(1), "H")
        self.assertEqual(constants.number_to_symbol(10), "Ne")

        with self.assertRaises(InvalidValueError):
            self.assertEqual(constants.number_to_symbol(0))

        with self.assertRaises(InvalidValueError):
            self.assertEqual(constants.number_to_symbol(1000))

    def test_symbol_to_mass(self):
        self.assertEqual(constants.symbol_to_mass("H"), 1.008)
        self.assertEqual(constants.symbol_to_mass("He"), 4.0026)
        self.assertEqual(constants.symbol_to_mass("O"), 15.999)
        self.assertEqual(constants.symbol_to_mass("P"), 30.974)
        self.assertEqual(constants.symbol_to_mass("H"), 1.008)
        self.assertEqual(constants.symbol_to_mass("Ne"), 20.180)

    def test_symbol_to_radius(self):
        self.assertEqual(constants.symbol_to_radius("H"), 0.53)
        self.assertEqual(constants.symbol_to_radius("He"), 0.31)
        self.assertEqual(constants.symbol_to_radius("O"), 0.48)
        self.assertEqual(constants.symbol_to_radius("P"), 0.98)
        self.assertEqual(constants.symbol_to_radius("H"), 0.53)
        self.assertEqual(constants.symbol_to_radius("Ne"), 0.38)

    def test_symbol_to_covalent_radius(self):
        self.assertEqual(constants.symbol_to_covalent_radius("H"), 0.37)
        self.assertEqual(constants.symbol_to_covalent_radius("He"), 0.32)
        self.assertEqual(constants.symbol_to_covalent_radius("O"), 0.73)
        self.assertEqual(constants.symbol_to_covalent_radius("P"), 1.06)
        self.assertEqual(constants.symbol_to_covalent_radius("H"), 0.37)
        self.assertEqual(constants.symbol_to_covalent_radius("Ne"), 0.69)

    def test_symbol_to_vdw_radius(self):
        self.assertEqual(constants.symbol_to_vdw_radius("H"), 1.20)
        self.assertEqual(constants.symbol_to_vdw_radius("He"), 1.40)
        self.assertEqual(constants.symbol_to_vdw_radius("O"), 1.52)
        self.assertEqual(constants.symbol_to_vdw_radius("P"), 1.80)
        self.assertEqual(constants.symbol_to_vdw_radius("H"), 1.20)
        self.assertEqual(constants.symbol_to_vdw_radius("Ne"), 1.54)

        with self.assertRaises(InvalidValueError):
            constants.symbol_to_vdw_radius("Sc")

    def test_symbol_to_free_polarizability(self):
        self.assertEqual(constants.symbol_to_free_polarizability("H"), constants.bohr_to_ang**3 * 4.50711)
        self.assertEqual(constants.symbol_to_free_polarizability("He"), constants.bohr_to_ang**3 * 1.38375)
        self.assertEqual(constants.symbol_to_free_polarizability("O"), constants.bohr_to_ang**3 * 5.3)
        self.assertEqual(constants.symbol_to_free_polarizability("P"), constants.bohr_to_ang**3 * 25)
        self.assertEqual(constants.symbol_to_free_polarizability("H"), constants.bohr_to_ang**3 * 4.50711)
        self.assertEqual(constants.symbol_to_free_polarizability("Ne"), constants.bohr_to_ang**3 * 2.66110)

    def test_symbol_to_ccsdt_free_polarizability(self):
        self.assertEqual(constants.symbol_to_ccsdt_free_polarizability("H"), 0.66582)
        self.assertEqual(constants.symbol_to_ccsdt_free_polarizability("B"), 3.64084)
        self.assertEqual(constants.symbol_to_ccsdt_free_polarizability("C"), 1.84613)
        self.assertEqual(constants.symbol_to_ccsdt_free_polarizability("N"), 1.08053)
        self.assertEqual(constants.symbol_to_ccsdt_free_polarizability("O"), 0.86381)
        self.assertEqual(constants.symbol_to_ccsdt_free_polarizability("H"), 0.66582)
        self.assertEqual(constants.symbol_to_ccsdt_free_polarizability("P"), 3.72507)

        with self.assertRaises(InvalidValueError):
            constants.symbol_to_ccsdt_free_polarizability("D")

        with self.assertRaises(InvalidValueError):
            constants.symbol_to_ccsdt_free_polarizability("Br")

suite = unittest.TestLoader().loadTestsFromTestCase(TestConstants)