import os
from potential_fitting.utils import files, SettingsReader
from potential_fitting.polynomials import MoleculeSymmetryParser
from .generate_mbnrg_fitting_code import generate_mbnrg_fitting_code

def prepare_mbnrg_fitting_code(settings_path, config_path, in_path, poly_path, poly_order, fit_path, use_direct):

    molecule = ""
    nmons = 0
    with open(in_path) as in_poly:
        while(True):
            fragment = in_poly.readline()
            if not fragment.startswith("add_molecule"):
                break

            if molecule is not "":
                molecule += "_"
            molecule += fragment[fragment.index("[") + 2:fragment.index("]") - 1]
            nmons += 1

    # copy needed files from poly_path to fit_path
    os.system("cp " + in_path + " " + fit_path + "/")

    # Ensure that if use_direct is active, the grad direct file exists
    directcpp_exists = os.path.isfile(poly_path + "/poly-direct.cpp") and os.path.isfile(poly_path + "/poly-grd-direct.cpp")

    if  use_direct and  not  directcpp_exists :
        raise FileDoesNotExistError("poly-(grd-)direct.cpp")
 
    degree = 0

    # find the number of variables and number of polynomials from the log file
    with open(poly_path + "/poly.log", "r") as poly_log:

        # loop thru each line
        for line in poly_log:

            # check if this line has the number of variables
            if "<> variables" in line:

                # parse number of variables from line
                number_of_variables = str(line[line.index("(") + 1: line.index(")")])
                break

        # loop thru each line
        for line in poly_log:

            # check if this line has a degree
            if "degree monomials" in line:

                # parse degree from line and update to max of all degrees so far
                # to find max degree of this polynomial
                degree = max(degree, int(line.split()[2]))

            # check if this line has the number of terms
            if "Total number of terms: " in line:

                # parse number of terms from line
                number_of_terms = str(line[line.index(":") + 2:])
                break

    config_obj = SettingsReader(config_path)
    config_obj.set("fitting","nvars",number_of_variables)
    config_obj.set("fitting","npoly",number_of_terms )
    config_obj.write(config_path)

    print("Executing python generator script")

    generate_mbnrg_fitting_code(settings_path, config_path, poly_order, in_path, poly_path, use_direct)

#    # restore settings
#    os.system("mv " + config_path + ".tmp " + config_path)

    # copy the template files
    os.system("cp " + os.path.dirname(os.path.abspath(__file__)) + "/fitting_code_templates/* " + fit_path)

    # move files from cwd into fit directory
    if (nmons<3):
        os.system("mv dispersion.* " + fit_path + "/")
        os.system("mv buckingham.* " + fit_path + "/")
    os.system("mv eval*b.cpp " + fit_path + "/")
    os.system("mv fit*b.cpp " + fit_path + "/")

    os.system("mv Makefile " + fit_path + "/")
    os.system("mv mbnrg_*_fit.* " + fit_path + "/")
    os.system("mv mon*.cpp mon*.h " + fit_path + "/")
    os.system("mv poly_*_fit.* " + fit_path + "/")

    # move the files for MBX into the directory MBX_files
    files.init_directory("MBX_files")
    symmetry_parser = MoleculeSymmetryParser("_".join(SettingsReader(settings_path).get("molecule", "symmetry").split(",")))
    system_name = symmetry_parser.get_symmetry()
    system_name = system_name.replace("(", "_o_").replace(")", "_c_")

    file_name = "mbnrg_{}b_{}_deg{}_v1.h".format(nmons, system_name, degree)
    os.system("mv " + file_name + " MBX_files/")
    file_name = "mbnrg_{}b_{}_deg{}_v1.cpp".format(nmons, system_name, degree)
    os.system("mv " + file_name + " MBX_files/")
    file_name = "poly_{}b_{}_deg{}_v1.h".format(nmons, system_name, degree)
    os.system("mv " + file_name + " MBX_files/")
    file_name = "poly_{}b_{}_deg{}_grad_v1.cpp".format(nmons, system_name, degree)
    os.system("mv " + file_name + " MBX_files/")
    file_name = "poly_{}b_{}_deg{}_nograd_v1.cpp".format(nmons, system_name, degree)
    os.system("mv " + file_name + " MBX_files/")
