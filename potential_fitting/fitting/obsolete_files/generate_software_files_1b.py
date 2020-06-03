import os, sys
import json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../../")
from potential_fitting.utils import SettingsReader

def write_cpp_software_properties(settings, config_file, fit_path, fit_cdl, mon_name, poly_order, molecule_in):
    # NOTE mon_name is the actual monomer name (h2o, ch4, c2h6...)

    # read config
    config = SettingsReader(config_file)

    # extract all the monomer information
    # Number of atoms in monomer
    nat = config.getint("fitting", "number_of_atoms")
    
    # Number of sites (electrostatic sites, NO LONE PAIRS)
    nsites = config.getint("fitting", "number_of_electrostatic_sites")
    
    # Obtain the lists with the excluded pairs
    excl12 = config.getlist("fitting", "excluded_pairs_12", int)[0]
    excl13 = config.getlist("fitting", "excluded_pairs_13", int)[0]
    excl14 = config.getlist("fitting", "excluded_pairs_14", int)[0]
    
    
    # Obtain charges (in the order of input), pols and polfacs
    chg = config.getlist("fitting", "charges", float)[0]
    pol = config.getlist("fitting", "polarizabilities", float)[0]
    polfac = config.getlist("fitting", "polarizability_fractions", float)[0]

    # Obtain polynomial coefficients and non_linear parameters
    constants = []
    polycoef = []
    npoly = -1
    print(fit_path + "/" + fit_cdl)
    pwd = os.getcwd()
    fit = open(pwd + "/" + fit_path + "/" + fit_cdl,'r')
    line = fit.readline()
    while True:
        # Skip name tag
        if line.strip().startswith(":"):
            if not "name" in line:
                constants.append(line.replace(":", "  m_"))
        # Find number of polynomial linear coefficients
        if line.startswith("  poly"):
            npoly = int(line.strip().split()[2].replace(";",""))
        # Store polynomial coefficients
        if line.startswith("poly"):
            for i in range(npoly):
                polycoef.append("            " + fit.readline().replace(";","};"))
        line = fit.readline()
        if line == "":
            break

    fit.close()

    # Open file that contains the code to add
    cppout = open("software_code.txt",'w')

    # Write code that needs to be added in the constructor
    cppout.write("=====>> SECTION CONSTRUCTOR <<=====\n")
    cppout.write("File: newly generated x1b...degN_v1x.cpp\n")
    a = '''
    if (mon == "''' + mon_name + '''") {
        coefficients = std::vector<double> { 
'''
    cppout.write(a)
    for c in polycoef:
        cppout.write(c)

    cppout.write("\n")
    
    for p in constants:
        cppout.write(p)

    a = '''
    } // end if mon == "''' + mon_name + '''" 


'''

    cppout.write(a)
    
    # Write code that needs to be added in the SITES section of the code
    cppout.write("=====>> SECTION SITES <<=====\n")
    cppout.write("File: src/bblock/sys_tools.cpp\n")
    a = '''
        } else if (mon[i] == "''' + mon_name + '''") {
            sites.push_back(''' + str(nsites) + ''');
            nat.push_back(''' + str(nat) + ''');



'''
    
    cppout.write(a)

    # Write code that needs to be added in the CHARGES section of the code
    cppout.write("=====>> SECTION CHARGES <<=====\n")
    cppout.write("File: src/bblock/sys_tools.cpp\n")

    a = '''
    } else if (mon_id == "''' + mon_name + '''") {
        for (size_t nv = 0; nv < n_mon; nv++) {
'''
    cppout.write(a)

    for i in range(len(chg)):
        cppout.write("            charges[fst_ind + nv*nsites + " + str(i) + "] = " + str(chg[i]) + " * CHARGECON;\n")

    a = '''
        }



'''

    cppout.write(a)

    # Write code that needs to be added in the POLFACS section of the code
    cppout.write("=====>> SECTION POLFACS <<=====\n")
    cppout.write("File: src/bblock/sys_tools.cpp\n")

    a = '''
    } else if (mon_id == "''' + mon_name + '''") {
        for (size_t nv = 0; nv < n_mon; nv++) {
'''
    cppout.write(a)

    for i in range(len(polfac)):
        cppout.write("            polfac[fst_ind + nv*nsites + " + str(i) + "] = " + str(polfac[i]) + ";\n")

    a = '''
        }



'''

    cppout.write(a)

    # Write code that needs to be added in the POLS section of the code
    cppout.write("=====>> SECTION POLS <<=====\n")
    cppout.write("File: src/bblock/sys_tools.cpp\n")

    a = '''
    } else if (mon_id == "''' + mon_name + '''") {
        for (size_t nv = 0; nv < n_mon; nv++) {
'''
    cppout.write(a)

    for i in range(len(pol)):
        cppout.write("            pol[fst_ind + nv*nsites + " + str(i) + "] = " + str(pol[i]) + ";\n")

    a = '''
        }



'''

    cppout.write(a)

    # Write code that needs to be added in excluded pair section
    cppout.write("=====>> SECTION EXCLUDED <<=====\n")
    cppout.write("File: src/bblock/sys_tools.cpp\n")

    a = '''
    if (mon == "''' + mon_name + '''") {
        // 12 distances
'''

    cppout.write(a)

    for i in excl12:
        cppout.write("        exc12.insert(std::make_pair(" + str(i[0]) + "," + str(i[1]) + "));\n")

    cppout.write("        // 13 distances\n")

    for i in excl13:
        cppout.write("        exc13.insert(std::make_pair(" + str(i[0]) + "," + str(i[1]) + "));\n")

    cppout.write("        // 14 distances\n")

    for i in excl14:
        cppout.write("        exc14.insert(std::make_pair(" + str(i[0]) + "," + str(i[1]) + "));\n")

    cppout.write("    }\n")

    # Write code that needs to be added in the ONEBODY_NOGRD section of the code
    cppout.write("=====>> SECTION ONEBODY_NOGRD <<=====\n")
    cppout.write("File: src/potential/1b/energy1b.cpp\n")

    a = '''
    } else if (mon == "''' + mon_name + '''") {
        x1b_''' + molecule_in + '''_deg''' + str(poly_order) + '''::x1b_''' + molecule_in + '''_v1x pot(mon);
        energies = pot.eval(xyz.data(), nm);


'''

    cppout.write(a)

    # Write code that needs to be added in the ONEBODY_GRD section of the code
    cppout.write("=====>> SECTION ONEBODY_GRD <<=====\n")
    cppout.write("File: src/potential/1b/energy1b.cpp\n")

    a = '''
    } else if (mon == "''' + mon_name + '''") {
        x1b_''' + molecule_in + '''_deg''' + str(poly_order) + '''::x1b_''' + molecule_in + '''_v1x pot(mon);
        energies = pot.eval(xyz.data(), grad.data(), nm, virial);


'''

    cppout.write(a)


    # Write code that needs to be added in the INCLUDE1B section of the code
    cppout.write("=====>> SECTION INCLUDE1B <<=====\n")
    cppout.write("File: src/potential/1b/energy1b.h\n")

    cppout.write("#include \"potential/1b/x1b_" + molecule_in + "_deg" + str(poly_order) + "_v1x.h\"\n")

    cppout.close()

    # software folder name
    sofdir = fit_path + "/" + "software_files"

    os.system("mv software_code.txt " + sofdir)


def generate_software_files_1b(settings, molecule_in, poly_path, poly_order, fit_path, config_file, fit_cdl, mon_name):
    """
    Generates the parts of the C++ code needed to add the PEF to the energy software

    Args:
        settings     - the file containing all relevent settings information
        molecule_in  - the A3B2 type molecule
        poly_path    - directory where polynomial files are
        poly_order   - the order of the polynomial in poly_path
        fit_path     - directory to generate fit code in
        config_file  - monomer config file
        fit_cdl      - the output cdl file of the fit to use
        mon_name     - the human understandable name of the monomer (co2, h2o, no3...)

    Returns:
        None
    """

    # software folder name
    sofdir = fit_path + "/" + "software_files"

    # create folder with the files for the software
    os.system("mkdir -p " + sofdir)

    # mv x1b files needed
    os.system("mv x1b_" + molecule_in + "_deg" + str(poly_order) +"_v1x.h " + sofdir)
    os.system("mv x1b_" + molecule_in + "_deg" + str(poly_order) +"_v1x.cpp " + sofdir)

    # define names for files
    headerf = sofdir + "/poly_1b_" + molecule_in + "_deg" + str(poly_order) + "_v1x.h"
    cppgrad = sofdir + "/poly_1b_" + molecule_in + "_deg" + str(poly_order) + "_v1x.cpp"
    cppnograd = sofdir + "/poly_1b_" + molecule_in + "_deg" + str(poly_order) + "_v1.cpp"

    # Copy polynomial files to temporary files
    os.system("cp " + poly_path + "/poly-model.h " + headerf + ".tmp")
    os.system("cp " + poly_path + "/poly-grd.cpp " + cppgrad + ".tmp")
    os.system("cp " + poly_path + "/poly-nogrd.cpp " + cppnograd + ".tmp")

    # modify header
    with open(headerf + ".tmp", "r") as grd:
        lines = grd.readlines()

    grd = open(headerf, "w")
    line_num = 0
    while line_num < len(lines):
        # go line by line, and replace keywords
        newline = lines[line_num].replace("POLY_MODEL","POLY_1B_" + molecule_in + "_DEG" + str(poly_order)).replace("mb_system","x1b_" + molecule_in + "_deg" + str(poly_order)).replace("poly-model","poly_1b_" + molecule_in + "_deg" + str(poly_order) + "_v1x").replace("poly_model","poly_1b_" + molecule_in + "_deg" + str(poly_order) + "_v1x")

        # skip eval_direct
        if not "eval_direct" in newline:
            grd.write(newline)
            line_num += 1
        else:
            line_num += 2
    
    # modify cpp grad file
    with open(cppgrad + ".tmp", "r") as grd:
        lines = grd.readlines()

    grd = open(cppgrad, "w")
    line_num = 0
    while line_num < len(lines):
        # go line by line, and replace keywords
        newline = lines[line_num].replace("poly_1b_" + molecule_in + "_v1x.h","poly_1b_" + molecule_in + "_deg" + str(poly_order)).replace("mb_system","x1b_" + molecule_in + "_deg" + str(poly_order)).replace("poly-model","poly_1b_" + molecule_in + "_deg" + str(poly_order) + "_v1x").replace("poly_model","poly_1b_" + molecule_in + "_deg" + str(poly_order) + "_v1x")
        
        line_num += 1
        grd.write(newline)

    # modify cpp nograd file
    with open(cppnograd + ".tmp", "r") as grd:
        lines = grd.readlines()

    grd = open(cppnograd, "w")
    line_num = 0
    while line_num < len(lines):
        # go line by line, and replace keywords
        newline = lines[line_num].replace("poly_1b_" + molecule_in + "_v1x.h","poly_1b_" + molecule_in + "_deg" + str(poly_order)).replace("mb_system","x1b_" + molecule_in + "_deg" + str(poly_order)).replace("poly-model","poly_1b_" + molecule_in + "_deg" + str(poly_order) + "_v1x").replace("poly_model","poly_1b_" + molecule_in + "_deg" + str(poly_order) + "_v1x")

        line_num += 1
        grd.write(newline)

    grd.close()

    # remove tmp files
    os.system("rm " + cppgrad + ".tmp " + cppnograd + ".tmp " + headerf + ".tmp")

    # Call the function to write the individual pieces of code, defined above
    write_cpp_software_properties(settings, config_file, fit_path, fit_cdl, mon_name, poly_order, molecule_in)
