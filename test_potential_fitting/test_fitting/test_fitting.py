import unittest
import os

from potential_fitting.utils import system, files
from potential_fitting import generate_mbnrg_fitting_code, compile_fit_code, prepare_fits, execute_fits, retrieve_best_fit

class TestFitting(unittest.TestCase):
    def setUpClass():

        TestFitting.resources_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
        TestFitting.out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

    def test_A2B1(self):

        molecule = "A1B2"
        config_path = os.path.join(TestFitting.resources_path, molecule, "config.ini")
        in_path = os.path.join(TestFitting.resources_path, molecule, "poly.in")
        poly_path = os.path.join(TestFitting.resources_path, molecule, "poly")
        settings_path = os.path.join(TestFitting.resources_path, molecule, "settings.ini")
        fitting_path = os.path.join(TestFitting.out_path, molecule, "fitting_code")
        fits_path = os.path.join(TestFitting.out_path, molecule, "fits")

        files.init_directory(fit_path)

        generate_mbnrg_fitting_code(settings_path, config_path, in_path, poly_path, 2, fitting_path, False)

        compile_fit_code(settings_path, fitting_path)

        training_path = os.path.join(TestFitting.resources_path, molecule, "training_set.xyz")
        fit_nc = "mbnrg.nc"
        fit_nc_path = os.path.join(TestFitting.out_path, molecule, "fits", fit_nc)

        prepare_fits(settings_path, fitting_path, 
                     training_path, fits_path, 
                     DE=20, alpha=0.0005, num_fits=1, 
                     ttm=False, over_ttm=False)
        execute_fits(settings_path, fits_path)
        retrieve_best_fit(settings_path, fits_path, fitted_nc_path = fit_nc)

        self.assertTrue(os.path.isfile(fit_nc_path))

#FIXME Add test A1
#    def test_A2B1_A2B1(self):
#
#        molecule = "A2B1_A2B1"
#        settings_path = os.path.join(TestFitting.resources_path, molecule, "settings.ini")
#        config_path = os.path.join(TestFitting.resources_path, molecule, "config.ini")
#        in_path = os.path.join(TestFitting.resources_path, molecule, "poly.in")
#        poly_path = os.path.join(TestFitting.resources_path, molecule, "poly")
#        ttm_path = os.path.join(TestFitting.out_path, molecule, "ttm")
#        fit_path = os.path.join(TestFitting.out_path, molecule, "fit")
#
#        files.init_directory(ttm_path)
#
#        generate_2b_ttm_fit_code(settings_path, config_path, molecule, ttm_path)
#
#        compile_fit_code(settings_path, ttm_path)
#
#        training_path = os.path.join(TestFitting.resources_path, molecule, "training_set.xyz")
#        ttm_fit_code_path = os.path.join(ttm_path, "fit-2b-ttm")
#
#        config_out_path = os.path.join(TestFitting.out_path, molecule, "config.ini")
#
#        system.call("cp", config_path, config_out_path)
#
#        fit_2b_ttm_training_set(settings_path, ttm_fit_code_path, training_path, ttm_path, config_out_path, 1)
#
#        files.init_directory(fit_path)
#        
#        prepare_2b_fitting_code(settings_path,
#                                config_out_path,
#                                in_path,
#                                poly_path,
#                                1,
#                                fit_path)
#
#
#        compile_fit_code(settings_path, fit_path)
#
#        fit_code_path = os.path.join(fit_path, "fit-2b")
#        fit_nc = os.path.join(TestFitting.out_path, molecule, "fit.nc")
#
#        fit_2b_training_set(settings_path, fit_code_path, training_path, fit_path, fit_nc, 1) 
#
#        self.assertTrue(os.path.isfile(fit_nc))





suite = unittest.TestLoader().loadTestsFromTestCase(TestFitting)
