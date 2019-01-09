import unittest

from potential_fitting.calculator import Model

class TestModel(unittest.TestCase):
	
	def test_get_method(self):
		model = Model("HF", "STO-3G", True)

		self.assertEqaul(model.get_method(), "HF")

		model = Model("wb97-mv", "cc-pvdz", False)

		self.assertEqaul(model.get_method(), "wb97-mv")

		model = Model("Nagic", "Wand", True)

		self.assertEqaul(model.get_method(), "Magic")

	def test_get_basis(self):
		model = Model("HF", "STO-3G", True)

		self.assertEqaul(model.get_basis(), "STO-3G")

		model = Model("wb97-mv", "cc-pvdz", False)

		self.assertEqaul(model.get_basis(), "cc-pvdz")

		model = Model("Magic", "Wand", True)

		self.assertEqaul(model.get_basis(), "Wand")

	def test_get_cp(self):
		model = Model("HF", "STO-3G", True)

		self.assertEqaul(model.get_cp(), True)

		model = Model("wb97-mv", "cc-pvdz", False)

		self.assertEqaul(model.get_cp(), False)

		model = Model("Magic", "Wand", True)

		self.assertEqaul(model.get_cp(), True)

	def test_equ(self):
		model1 = Model("HF", "STO-3G", True)
		
		model2 = Model("HF", "STO-3G", True)

		self.assertTrue(model1 == model2)

		model1 = Model("hf", "sto-3g", True)

		model2 = Model("HF", "STO-3G", True)

		self.assertFalse(model1 == model2)	

suite = unittest.TestLoader().loadTestsFromTestCase(TestAtom)