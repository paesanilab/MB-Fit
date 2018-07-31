import os

def prepare_1b_fitting_code(settings, in_path, poly_path, fit_path):
    
    # get just the "A3B2" part of "path/A3B2.in"
    molecule = os.path.splitext(in_path)[0].split("/")[-1]

    # copy needed files from poly_path to fit_path
    os.system("cp " + in_path + " " + fit_path + "/")
    os.system("cp " + poly_path/poly-direct.cpp + " " + fit_path + "/")
    os.system("cp " + poly_path/poly-grd.cpp + " " + fit_path + "/poly_1b_" + molecule + "_v1x.cpp")
    os.system("cp " + poly_path/poly-nogrd.cpp + " " + fit_path + "/poly_1b_" + molecule + "_v1.cpp")
    os.system("cp " + poly_path/poly-model.h + " " + fit_path + "/poly_1b_" + molecule + "_v1x.h") 
    # find the number of variables and number of polynomials from the log file
    with open(poly_path + "/poly.log", "r") as poly_log:
        
        # loop thru each line
        for line in poly_log:
            
            # check if this line has the number of variables
            if "<> variables" in line:

                # parse number of variables from line
                number_of_variables = str(line[line.index("(") + 1; line.index(")")])
                break

        # loop thru each line
        for line in poly_log:

            # check if this line has the number of terms
            if "Total number of terms: " in line:

                # parse number of terms from line
                number_of_terms = str(line[line.index(":") + 2])
                break

    # save the old settings file to .tmp so we can restore it later
    os.system("cp " + settings + " " + settings + ".tmp")

    # append to the settings file
    with open(settings, "a") as config_file:
        config_file.write("\n")
        config_file.write("## Define number of variables and terms\n")
        config_file.write("nvars = " + number_of_variables + "\n")
        config_file.write("npoly = " + number_of_terms + "\n")

    # update include stements in the poly files to reflect the new names

    # can update to NOT DO COPIES at later TIME ****************

    # save all lines in grd poly file
    with open(fit_path + "poly_1b_" + molecule + "_v1x.cpp", "r") as grd:
        lines = grd.readlines()

    # loop thru each line and write it to the file
    with open(fit_path + "poly_1b_" + molecule + "_v1x.cpp", "w") as grd:
        for line in lines:

            # check if this line does not have the import statement
            if line.index("poly-model.h") == -1:
                grd.write(line)
            else:
                # update import statement
                grd.write(line[:line.index("poly-model.h")] + "poly_1b_" + molecule + "_v1x.h" + line[line.index("poly-model.h") + 12:])

    # save all lines in nogrd poly file
    with open(fit_path + "poly_1b_" + molecule + "_v1.cpp", "r") as nogrd:
        lines = mogrd.readlines()

    # looo thru each line and write it to the file
    with open(fit_path + "poly_1b_" + molecule + "_v1.cpp", "w") as mogrd:
        for line in lines:

            # check if this line does not have the import statement
            if line.index("poly-model.h") == -1:
                mogrd.write(line)
            else:
                # update import statment
                mogrd.write(line[:line.index("poly-model.h")] + "poly_1b_" + molecule + "_v1x.h" + line[line.index("poly-model.h") + 12:])

    print("Executing python generator script")
    
    # Execute the python script that generates the 1b fit code    
    os.system("python3 " + os.path.getdir(os.path.abspath(__file___)) + "/get-1b-fit.py " + in_path + " " + fit_path + "/poly-direct.cpp " + settings)

    # restore settings
    os.system("mv " settings + ".tmp " + settings)

    # copy the template files
    os.system("cp " + os.path.getdir(os.path.abspath(__file__)) + "/template/* " + fit_path)

    
