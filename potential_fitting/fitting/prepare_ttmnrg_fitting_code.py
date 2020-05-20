import os
from potential_fitting.utils import files, SettingsReader
from potential_fitting.polynomials import MoleculeSymmetryParser
from .generate_ttmnrg_fitting_code import generate_ttmnrg_fitting_code

def prepare_ttmnrg_fitting_code(settings_path, config_path, fit_path):

    molecule = ""
    nmons = 0

    config_obj = SettingsReader(config_path)

    print("Executing python generator script")

    generate_ttmnrg_fitting_code(settings_path, config_path)

    # copy the template files
    os.system("cp " + os.path.dirname(os.path.abspath(__file__)) + "/fitting_code_templates/* " + fit_path)

    # move files from cwd into fit directory
    os.system("mv dispersion.* " + fit_path + "/")
    os.system("mv buckingham.* " + fit_path + "/")
    os.system("mv eval*b.cpp " + fit_path + "/")
    os.system("mv fit*b.cpp " + fit_path + "/")

    os.system("mv eval*b-ttm.cpp " + fit_path + "/")
    os.system("mv fit*b-ttm.cpp " + fit_path + "/")

    os.system("mv Makefile " + fit_path + "/")
    os.system("mv mon*.cpp mon*.h " + fit_path + "/")
