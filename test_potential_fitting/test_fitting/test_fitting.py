import unittest
import os

from potential_fitting.utils import system, files
from potential_fitting import compile_fit_code, fit_1b_training_set, fit_2b_ttm_training_set, fit_2b_training_set, generate_2b_ttm_fit_code
from potential_fitting.fitting import prepare_1b_fitting_code, prepare_2b_fitting_code

class TestFitting(unittest.TestCase):
    def setUpClass():

        TestFitting.resources_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
        TestFitting.out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

    def test_A2B1(self):

        molecule = "A2B1"
        config_path = os.path.join(TestFitting.resources_path, molecule, "config.ini")
        in_path = os.path.join(TestFitting.resources_path, molecule, "poly.in")
        poly_path = os.path.join(TestFitting.resources_path, molecule, "poly")
        fit_path = os.path.join(TestFitting.out_path, molecule, "fit")

        files.init_directory(fit_path)

        prepare_1b_fitting_code(config_path,
                                molecule,
                                in_path,
                                poly_path,
                                2,
                                fit_path)

        settings_path = os.path.join(TestFitting.resources_path, molecule, "settings.ini")

        compile_fit_code(settings_path, fit_path)

        training_path = os.path.join(TestFitting.resources_path, molecule, "training_set.xyz")
        fit_code_path = os.path.join(fit_path, "fit-1b")
        fit_nc = os.path.join(TestFitting.out_path, molecule, "fit.nc")

        fit_1b_training_set(settings_path, fit_code_path, training_path, fit_path, fit_nc, 1) 

        self.assertTrue(os.path.isfile(fit_nc))

    def test_A2B1_A2B1(self):

        molecule = "A2B1_A2B1"
        settings_path = os.path.join(TestFitting.resources_path, molecule, "settings.ini")
        config_path = os.path.join(TestFitting.resources_path, molecule, "config.ini")
        in_path = os.path.join(TestFitting.resources_path, molecule, "poly.in")
        poly_path = os.path.join(TestFitting.resources_path, molecule, "poly")
        ttm_path = os.path.join(TestFitting.out_path, molecule, "ttm")
        fit_path = os.path.join(TestFitting.out_path, molecule, "fit")

        files.init_directory(ttm_path)

        generate_2b_ttm_fit_code(settings_path, config_path, molecule, ttm_path)

        compile_fit_code(settings_path, ttm_path)

        training_path = os.path.join(TestFitting.resources_path, molecule, "training_set.xyz")
        ttm_fit_code_path = os.path.join(ttm_path, "fit-2b-ttm")

        config_out_path = os.path.join(TestFitting.out_path, molecule, "config.ini")

        system.call("cp", config_path, config_out_path)

        fit_2b_ttm_training_set(settings_path, ttm_fit_code_path, training_path, ttm_path, config_out_path, 1)

        files.init_directory(fit_path)
        
        prepare_2b_fitting_code(settings_path,
                                config_out_path,
                                in_path,
                                poly_path,
                                1,
                                fit_path)


        compile_fit_code(settings_path, fit_path)

        fit_code_path = os.path.join(fit_path, "fit-2b")
        fit_nc = os.path.join(TestFitting.out_path, molecule, "fit.nc")

        fit_2b_training_set(settings_path, fit_code_path, training_path, fit_path, fit_nc, 1) 

        self.assertTrue(os.path.isfile(fit_nc))





suite = unittest.TestLoader().loadTestsFromTestCase(TestFitting)
