import unittest

from potential_fitting.calculator import Model

class TestModel(unittest.TestCase):

    # set up before the first test case
    def setUpClass():
        pass

    # clean up after the last test case
    def tearDownClass():
        pass

    # set up before each test case
    def setUp(self):

        self.HF_STO3G_True = Model("HF", "STO-3G", True)
        self.wb97mv_ccpvdz_False = Model("wb97m-v", "cc-pvdz", False)
        self.magic_wand = Model("Magic", "Wand", True)

    # clean up after each test case
    def tearDown(self):
        pass
    
    def test_get_method(self):

        self.assertEqaul(self.HF_STO3G_True.get_method(), "HF")
        self.assertEqaul(self.wb97mv_ccpvdz_false.get_method(), "wb97m-v")
        self.assertEqaul(self.magic_wand.get_method(), "Magic")

    def test_get_basis(self):
        
        self.assertEqaul(self.HF_STO3G_True.get_basis(), "HF")
        self.assertEqaul(self.wb97mv_ccpvdz_false.get_basis(), "cc-pvdz")
        self.assertEqaul(self.magic_wand.get_basis(), "Wand")

    def test_get_cp(self):
        
        self.assertEqaul(self.HF_STO3G_True.get_cp(), True)
        self.assertEqaul(self.wb97mv_ccpvdz_false.get_cp(), False)
        self.assertEqaul(self.magic_wand.get_cp(), True)

    def test_eq(self):

        model1 = Model("HF", "STO-3G", True)
        self.assertTrue(model1 == self.HF_STO3G_True)

        model2 = Model("hf", "sto-3g", True)
        self.assertFalse(model2 == self.HF_STO3G_True)

        model3 = Model("wb97", "m-vcc-pvdz", False)
        self.assertFalse(model3 == self.wb97mv_ccpvdz_False)

    def test_ne(self):

        model1 = Model("HF", "STO-3G", True)
        self.assertFalse(model1 != self.HF_STO3G_True)

        model2 = Model("hf", "sto-3g", True)
        self.assertTrue(model2 != self.HF_STO3G_True)

        model3 = Model("wb97", "m-vcc-pvdz", False)
        self.assertTrue(model3 != self.wb97mv_ccpvdz_False)

suite = unittest.TestLoader().loadTestsFromTestCase(TestModel)