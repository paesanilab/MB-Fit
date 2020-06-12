import unittest
import os

from test_potential_fitting.test_case_with_id import TestCaseWithId
from potential_fitting.utils import system, files
from potential_fitting import generate_mbnrg_fitting_code, compile_fit_code, prepare_fits, execute_fits, retrieve_best_fit, generate_ttmnrg_fitting_code, update_config_with_ttm

class TestFitting(TestCaseWithId):
    def setUpClass():

        TestFitting.resources_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "resources")
        TestFitting.out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

    def test_A1B2(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))

        molecule = "A1B2"
        config_in_path = os.path.join(TestFitting.resources_path, molecule, "config.ini")
        config_path = os.path.join(TestFitting.out_path, molecule, "config.ini")
        in_path = os.path.join(TestFitting.resources_path, molecule, "poly.in")
        poly_path = os.path.join(TestFitting.resources_path, molecule, "poly")
        settings_path = os.path.join(TestFitting.resources_path, molecule, "settings.ini")
        fitting_path = os.path.join(TestFitting.out_path, molecule, "fitting_code")
        fits_path = os.path.join(TestFitting.out_path, molecule, "fits")

        files.init_directory(os.path.join(TestFitting.out_path, molecule))
        os.system("cp " + config_in_path + " " + config_path)

        files.init_directory(fits_path)

        generate_mbnrg_fitting_code(settings_path, config_path, in_path, poly_path, 2, fitting_path, False)

        compile_fit_code(settings_path, fitting_path)

        training_path = os.path.join(TestFitting.resources_path, molecule, "training_set.xyz")
        fit_nc = "mbnrg.nc"
        fit_nc_path = os.path.join(fits_path, "best_fit", fit_nc)

        prepare_fits(settings_path, fitting_path, 
                     training_path, fits_path, 
                     DE=20, alpha=0.0005, num_fits=1, 
                     ttm=False, over_ttm=False)
        execute_fits(settings_path, fits_path)
        retrieve_best_fit(settings_path, fits_path, fitted_nc_path = fit_nc)
        
        self.assertTrue(os.path.isfile(fit_nc_path))

        self.test_passed = True

    def test_A1B2_A1B2_ttmnrg(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        molecule = "A1B2_A1B2"
        settings_path = os.path.join(TestFitting.resources_path, molecule + "_TTM", "settings.ini")
        config_in_path = os.path.join(TestFitting.resources_path, molecule + "_TTM", "config.ini")
        config_path = os.path.join(TestFitting.out_path, molecule + "_TTM", "config.ini")
        ttm_fit_code_path = os.path.join(TestFitting.out_path, molecule + "_TTM", "ttm")
        ttm_fits_path = os.path.join(TestFitting.out_path, molecule + "_TTM", "fit")
        training_set = os.path.join(TestFitting.resources_path, molecule + "_TTM","training_set.xyz")
        ttmparams = os.path.join(ttm_fits_path, "best_fit", "ttm-nrg_params.dat")

        files.init_directory(os.path.join(TestFitting.out_path, molecule + "_TTM"))
        os.system("cp " + config_in_path + " " + config_path)

        files.init_directory(ttm_fits_path)
        
        generate_ttmnrg_fitting_code(settings_path, config_path, ttm_fit_code_path)
        compile_fit_code(settings_path, ttm_fit_code_path)

        prepare_fits(settings_path, ttm_fit_code_path, 
                               training_set, ttm_fits_path, 
                               DE=20, alpha=0.0005, num_fits=2, 
                               ttm=True, over_ttm=False)

        execute_fits(settings_path, ttm_fits_path)

        retrieve_best_fit(settings_path, ttm_fits_path)

        update_config_with_ttm(settings_path, ttm_fits_path, config_path)

        self.assertTrue(os.path.isfile(ttmparams))

        self.test_passed = True

    def test_A1B2_A1B2_mbnrg(self):
        self.test_folder = os.path.dirname(os.path.abspath(__file__))
        molecule = "A1B2_A1B2"
        config_in_path = os.path.join(TestFitting.resources_path, molecule + "_MB", "config.ini")
        config_path = os.path.join(TestFitting.out_path, molecule + "_MB", "config.ini")
        in_path = os.path.join(TestFitting.resources_path, molecule + "_MB", "poly.in")
        poly_path = os.path.join(TestFitting.resources_path, molecule + "_MB", "poly")
        settings_path = os.path.join(TestFitting.resources_path, molecule + "_MB", "settings.ini")
        fitting_path = os.path.join(TestFitting.out_path, molecule + "_MB", "fitting_code")
        fits_path = os.path.join(TestFitting.out_path, molecule + "_MB", "fits")

        files.init_directory(os.path.join(TestFitting.out_path, molecule + "_MB"))
        os.system("cp " + config_in_path + " " + config_path)

        files.init_directory(fits_path)

        generate_mbnrg_fitting_code(settings_path, config_path, in_path, poly_path, 2, fitting_path, False)

        compile_fit_code(settings_path, fitting_path)

        training_path = os.path.join(TestFitting.resources_path, molecule + "_MB", "training_set.xyz")
        fit_nc = "mbnrg.nc"
        fit_nc_path = os.path.join(TestFitting.out_path, molecule + "_MB", "fits", "best_fit", fit_nc)

        prepare_fits(settings_path, fitting_path,
                     training_path, fits_path,
                     DE=20, alpha=0.0005, num_fits=1,
                     ttm=False, over_ttm=False)
        execute_fits(settings_path, fits_path)
        retrieve_best_fit(settings_path, fits_path, fitted_nc_path = fit_nc)

        self.assertTrue(os.path.isfile(fit_nc_path))

        self.test_passed = True
        
suite = unittest.TestLoader().loadTestsFromTestCase(TestFitting)
