import os

def prepare_fitting_code(settings_path, config_path, in_path, poly_path, poly_order, fit_path):
    
    # get just the "A3B2" part of "path/A3B2.in"
    molecule = ""
    with open(in_path) as in_poly:
        while(True):
            fragment = in_poly.readline()
            if not fragment.startswith("add_molecule"):
                break
        
            if molecule is not "":
                molecule += "_"
            molecule += fragment[fragment.index("[") + 2:fragment.index("]") - 1]

    # copy needed files from poly_path to fit_path
    os.system("cp " + in_path + " " + fit_path + "/")
    os.system("cp " + poly_path + "/poly-direct.cpp " + fit_path + "/")
#    os.system("cp " + poly_path + "/poly-grd.cpp " + fit_path + "/poly_2b_" + molecule + "_v1x.cpp")
#    os.system("cp " + poly_path + "/poly-nogrd.cpp " + fit_path + "/poly_2b_" + molecule + "_v1.cpp")
#    os.system("cp " + poly_path + "/poly-model.h " + fit_path + "/poly_2b_" + molecule + "_v1x.h") 
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

            # check if this line has the number of terms
            if "Total number of terms: " in line:

                # parse number of terms from line
                number_of_terms = str(line[line.index(":") + 2:])
                break

    # save the old settings file to .tmp so we can restore it later
    os.system("cp " + config_path + " " + config_path + ".tmp")

    # append to the settings file
    with open(config_path, "a") as config_file:
        config_file.write("\n")
        config_file.write("## Define number of variables and terms\n")
        config_file.write("nvars = " + number_of_variables + "\n")
        config_file.write("npoly = " + number_of_terms + "\n")

    # update include stements in the poly files to reflect the new names

    # can update to NOT DO COPIES at later TIME ****************

    # save all lines in grd poly file
#    with open(fit_path + "/poly_2b_" + molecule + "_v1x.cpp", "r") as grd:
#        lines = grd.readlines()
#
#    # loop thru each line and write it to the file
#    with open(fit_path + "/poly_2b_" + molecule + "_v1x.cpp", "w") as grd:
#        for line in lines:
#
#            # check if this line has the import statement
#            try:
#                line.index("poly-model.h")
#                # update import statement
#                grd.write(line[:line.index("poly-model.h")] + "poly_2b_" + molecule + "_v1x.h" + line[line.index("poly-model.h") + 12:])
#            except ValueError:
#                # write original line if no import statement
#                grd.write(line)
#
#    # save all lines in nogrd poly file
#    with open(fit_path + "/poly_2b_" + molecule + "_v1.cpp", "r") as nogrd:
#        lines = nogrd.readlines()
#
#    # looo thru each line and write it to the file
#    with open(fit_path + "/poly_2b_" + molecule + "_v1.cpp", "w") as nogrd:
#        for line in lines:
#
#            # check if this line has the import statement
#            try:
#                line.index("poly-model.h")
#                # update import statment
#                nogrd.write(line[:line.index("poly-model.h")] + "poly_2b_" + molecule + "_v1x.h" + line[line.index("poly-model.h") + 12:])
#            except ValueError:
#                # write original line of no import statement
#                nogrd.write(line)

    print("Executing python generator script")
    
    # Execute the python script that generates the 1b fit code    
    os.system("python3 " + os.path.dirname(os.path.abspath(__file__)) + "/../../codes/nb-codes/generate_fitting_code.py " + settings_path + " " + config_path + " " + fit_path + "/poly-direct.cpp " + str(poly_order) + " " + in_path)

    # restore settings
    os.system("mv " + config_path + ".tmp " + config_path)

    # copy the template files
    os.system("cp " + os.path.dirname(os.path.abspath(__file__)) + "/../../codes/nb-codes/template/* " + fit_path)

    # move files from cwd into fit directory
    os.system("mv dispersion.* " + fit_path + "/")
    os.system("mv buckingham.* " + fit_path + "/")
    os.system("mv eval*b.cpp " + fit_path + "/")
    os.system("mv fit*b.cpp " + fit_path + "/")
    os.system("mv Makefile " + fit_path + "/")
    os.system("mv mbnrg_*_fit.* " + fit_path + "/")
    os.system("mv mon*.cpp mon*.h " + fit_path + "/")
    os.system("mv poly_*_fit.* " + fit_path + "/")

#    os.system("mv buckingham.h " + fit_path + "/")
#    os.system("mv fit-2b-wbuck.cpp " + fit_path + "/")
#    os.system("mv fit-2b.cpp " + fit_path + "/")
#    os.system("mv mon1.h " + fit_path + "/")
#    os.system("mv mon2.h " + fit_path + "/")
#    os.system("mv poly_2b_" + molecule +".h " + fit_path + "/")
#    os.system("mv x2b_" + molecule + "_v1.cpp " + fit_path + "/")
#    os.system("mv x2b_" + molecule + "_v1x.cpp " + fit_path + "/")
#    os.system("mv buckingham.cpp " + fit_path + "/")
#    os.system("mv eval-2b.cpp " + fit_path + "/")
#    os.system("mv eval-2b-wbuck.cpp " + fit_path + "/")
#    os.system("mv mon1.cpp " + fit_path + "/")
#    os.system("mv mon2.cpp " + fit_path + "/")
#    os.system("mv poly_2b_" + molecule + ".cpp " + fit_path + "/")
#    os.system("mv training_set.h " + fit_path + "/")
#    os.system("mv x2b_" + molecule + "_v1.h " + fit_path + "/")
#    os.system("mv x2b_" + molecule + "_v1x.h " + fit_path + "/")
