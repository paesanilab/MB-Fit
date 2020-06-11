import unittest, os

from potential_fitting.calculator import Model

class TestModel(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestModel, self).__init__(*args, **kwargs)
        self.test_passed = False
        self.test_name = self.id()

    # set up before each test case
    def setUp(self):

        self.HF_STO3G_True = Model("HF", "STO-3G", True)
        self.wb97mv_ccpvdz_False = Model("wb97m-v", "cc-pvdz", False)
        self.magic_wand = Model("Magic", "Wand", True)

    # clean up after each test case
    def tearDown(self):
        local_output = os.path.join(os.path.dirname(os.path.abspath(__file__)),"output")
        mbfithome = os.environ.get('MBFIT_HOME')
        if os.path.isdir(local_output):
            if self.test_passed:
                os.system("mkdir -p " + os.path.join(mbfithome, "passed_tests_outputs"))
                os.system("mv " + local_output + " " + os.path.join(mbfithome, "passed_tests_outputs", self.test_name))
            else:
                os.system("mkdir -p " + os.path.join(mbfithome, "failed_tests_outputs"))
                os.system("mv " + local_output + " " + os.path.join(mbfithome, "failed_tests_outputs", self.test_name))
    
    def test_get_method(self):

        self.assertEqual(self.HF_STO3G_True.get_method(), "HF")
        self.assertEqual(self.wb97mv_ccpvdz_False.get_method(), "wb97m-v")
        self.assertEqual(self.magic_wand.get_method(), "Magic")

        self.test_passed = True

    def test_get_basis(self):
        
        self.assertEqual(self.HF_STO3G_True.get_basis(), "STO-3G")
        self.assertEqual(self.wb97mv_ccpvdz_False.get_basis(), "cc-pvdz")
        self.assertEqual(self.magic_wand.get_basis(), "Wand")

        self.test_passed = True

    def test_get_cp(self):
        
        self.assertEqual(self.HF_STO3G_True.get_cp(), True)
        self.assertEqual(self.wb97mv_ccpvdz_False.get_cp(), False)
        self.assertEqual(self.magic_wand.get_cp(), True)

        self.test_passed = True

    def test_eq(self):

        model1 = Model("HF", "STO-3G", True)
        self.assertTrue(model1 == self.HF_STO3G_True)

        model2 = Model("hf", "sto-3g", True)
        self.assertFalse(model2 == self.HF_STO3G_True)

        model3 = Model("wb97", "m-vcc-pvdz", False)
        self.assertFalse(model3 == self.wb97mv_ccpvdz_False)

        self.test_passed = True

    def test_ne(self):

        model1 = Model("HF", "STO-3G", True)
        self.assertFalse(model1 != self.HF_STO3G_True)

        model2 = Model("hf", "sto-3g", True)
        self.assertTrue(model2 != self.HF_STO3G_True)

        model3 = Model("wb97", "m-vcc-pvdz", False)
        self.assertTrue(model3 != self.wb97mv_ccpvdz_False)

        self.test_passed = True

suite = unittest.TestLoader().loadTestsFromTestCase(TestModel)
