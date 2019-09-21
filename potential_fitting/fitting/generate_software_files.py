import os, sys
import json
import math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../../")
from potential_fitting.utils import SettingsReader
from potential_fitting.polynomials import FragmentParser
from potential_fitting.exceptions import InconsistentValueError

# FIXME duplicated function in utils_nb_fitting.
def get_atom_types(fragment):
    fragment_parser = FragmentParser(fragment, 'a')

    atom_list = []

    for atom_type in fragment_parser.get_atom_and_virtual_site_types():
        atom_list.append(atom_type.get_type())
        atom_list.append(atom_type.get_count())

    return atom_list

def generate_software_files(settings_path, config_file, mon_ids, degree, ttm_only = False, MBX_HOME = None, version = "v1"):
    # NOTE mon_ids is a list with actual monomer names (h2o, ch4, c2h6...)
    # Read Settings
    settings = SettingsReader(settings_path)
    # Read config
    config = SettingsReader(config_file)

    ############################################################################
    ## MONOMER PROPERTIES ######################################################
    ############################################################################

    # Obtain monomers from settings
    monomers = settings.get("molecule", "symmetry").split(",")
    number_of_monomers = len(monomers)
    
    # Store number of atoms in each monomer in nat1 and nat2
    number_of_atoms = [int(i) for i in settings.get("molecule", "fragments").split(",")]
    
    # Set the number of sites
    number_of_sites = number_of_atoms
    
    # Get which monomer should be mb-pol monomer (if any)
    use_mbpol = [int(i) for i in settings.get("molecule", "use_mbpol").split(",")]
    
    # Define if lone pairs are used based on monomer names
    # Update number of sites if needed
    use_lonepairs = [0]*number_of_monomers
    for i in range(number_of_monomers):
        if "X" in monomers[i] or "Y" in monomers[i] or "Z" in monomers[i]:
            use_lonepairs[i] = 1
        if use_mbpol[i] != 0:
            number_of_sites[i] += 1
    
    # Obtain the lists with the excluded pairs
    excluded_pairs_12 = config.getlist("fitting", "excluded_pairs_12", int)
    excluded_pairs_13 = config.getlist("fitting", "excluded_pairs_13", int)
    excluded_pairs_14 = config.getlist("fitting", "excluded_pairs_14", int)
    
    #Obtain charges (in the order of input), pols and polfacs
    charges = config.getlist("fitting", "charges", float)
    polarizabilities = config.getlist("fitting", "polarizabilities", float)
    polarizability_factors = config.getlist("fitting", "polarizability_factors", float)

    ############################################################################
    ## SYSTEM  PROPERTIES ######################################################
    ############################################################################

    # Obtain A, C6 and d6 from user in the same order as the given pairs AA, AB ...
    # By default, if not specified in the config.ini file, these lists
    # will be set to a list of the same length, but with all coefficients 0.0
    C6 = config.getlist("fitting", "c6", float)
    
    # last element of A_constants is list of the inter_molecular A constants
    # Initialize A and b to 0, for now
    try:
        A_buck = config.getlist("fitting", "A", float)
        b_buck = config.getlist("fitting", "d6", float)
        d6 = b_buck
    except:
        A_buck = [0.0] * len(C6)
        b_buck = [0.0] * len(C6)
        d6 = [0.0] * len(C6)
        # Throw warning if 1b or 2b. This should be defined by the time 
        # this is called in those cases

    ############################################################################
    ## Polynomial ##############################################################
    ############################################################################
    
    print("Getting polynomial fitted parameters...")
    # Define a couple things
    workdir = os.getcwd()
    mbnrg_best_fit = workdir + "/" + settings.get("files", "log_path") + "/mb_nrg_fits/best_fit"
    cdl_file = "fit-" + str(number_of_monomers) + "b.cdl"

    # Obtain polynomial coefficients and non_linear parameters
    constants = []
    polycoef = []
    npoly = -1
    
    cdl = open(mbnrg_best_fit + "/" + cdl_file,'r')
    line = cdl.readline()
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
                polycoef.append("            " + cdl.readline().replace(";","};"))
        line = cdl.readline()
        if line == "":
            break

    cdl.close()

    my_constructor_text = ""
    mon_id_sorted = sorted(list(enumerate(mon_ids)), key=lambda x: x[1])

    my_constructor_text += "    if ("
    ids = []
    for i in range(len(mon_id_sorted)):
        ids.append("mon" + str(i+1) + " == \"" + mon_ids[i] + "\"")

    my_constructor_text += " and ".join(ids) + ") {\n"

    ## Open file that contains the code to add
    #cppout = open("software_code.txt",'w')

    ## Write code that needs to be added in the constructor
    #cppout.write("=====>> SECTION CONSTRUCTOR <<=====\n")
    #cppout.write("File: newly generated x1b...degN_v1x.cpp\n")
    #a = '''
    for c in polycoef:
        my_constructor_text += c

    c += "\n"
    
    for p in constants:
        my_constructor_text += p

    my_constructor_text += "\n    } // end if " + " and ".join(ids) + "\n" 

    ############################################################################
    ## Monomer Properties ##### Only if one-body ###############################
    ############################################################################
    i = 100
    if number_of_monomers == 1:
        print("Getting monomer properties...")
    # Write code that needs to be added in the SITES section of the code
    #cppout.write("=====>> SECTION SITES <<=====\n")
    #cppout.write("File: src/bblock/sys_tools.cpp\n")

        my_monomer_sites_text = ""
        my_monomer_sites_text += """
        }} else if (mon[i] == "{}") {{
            sites.push_back({});
            nat.push_back({});
""".format(mon_ids[0],number_of_sites[0],number_of_atoms[0])

    
    # Write code that needs to be added in the CHARGES section of the code
    #cppout.write("=====>> SECTION CHARGES <<=====\n")
    #cppout.write("File: src/bblock/sys_tools.cpp\n")

        my_monomer_charges_text = ""
        my_monomer_charges_text += """
        }} else if (mon_id == "{}") {{
            for (size_t nv = 0; nv < n_mon; nv++) {{
""".format(mon_ids[0])

        for j in range(len(charges[0])):
            my_monomer_charges_text += "                charges[fst_ind + nv*nsites + {}] = {} * CHARGECON;\n".format(j,charges[0][j])
        my_monomer_charges_text += "            }\n"


    # Write code that needs to be added in the POLFACS section of the code
    #cppout.write("=====>> SECTION POLFACS <<=====\n")
    #cppout.write("File: src/bblock/sys_tools.cpp\n")

        my_monomer_polfacs_text = ""
        my_monomer_polfacs_text += """
        }} else if (mon_id == "{}") {{
            for (size_t nv = 0; nv < n_mon; nv++) {{
""".format(mon_ids[0])

        for j in range(len(polarizability_factors[0])):
            my_monomer_polfacs_text += "                polfac[fst_ind + nv*nsites + {}] = {};\n".format(j,polarizability_factors[0][j])
        my_monomer_polfacs_text += "            }\n"


    # Write code that needs to be added in the POLS section of the code
    #cppout.write("=====>> SECTION POLS <<=====\n")
    #cppout.write("File: src/bblock/sys_tools.cpp\n")

        my_monomer_pols_text = ""
        my_monomer_pols_text += """
        }} else if (mon_id == "{}") {{
            for (size_t nv = 0; nv < n_mon; nv++) {{
""".format(mon_ids[0])

        for j in range(len(polarizabilities[0])):
            my_monomer_pols_text += "                pol[fst_ind + nv*nsites + {}] = {};\n".format(j,polarizabilities[0][j])
        my_monomer_pols_text += "            }\n"


    # Write code that needs to be added in excluded pair section
    #cppout.write("=====>> SECTION EXCLUDED <<=====\n")
    #cppout.write("File: src/bblock/sys_tools.cpp\n")

        my_monomer_excluded_text = ""
        my_monomer_excluded_text += """
        if (mon == "{}") {{
            // 12 distances
""".format(mon_ids[0])

        for excl in excluded_pairs_12[0]:
            my_monomer_excluded_text += "            exc12.insert(std::make_pair({}, {}));\n".format(excl[0],excl[1])

        my_monomer_excluded_text += "            // 13 distances\n"

        for excl in excluded_pairs_13[0]:
            my_monomer_excluded_text += "            exc13.insert(std::make_pair({}, {}));\n".format(excl[0],excl[1])

        my_monomer_excluded_text += "            // 14 distances\n"

        for excl in excluded_pairs_14[0]:
            my_monomer_excluded_text += "            exc14.insert(std::make_pair({}, {}));\n".format(excl[0],excl[1])

        my_monomer_excluded_text += "        }\n"

    # Write code that needs to be added in the ONEBODY_NOGRD section of the code
    #cppout.write("=====>> SECTION ONEBODY_NOGRD <<=====\n")
    #cppout.write("File: src/potential/1b/energy1b.cpp\n")

    print("Getting energy calls...")
    system_name = monomers[0]
    for i in range(1,number_of_monomers):
        system_name += "_" + monomers[i]

    ids = []
    for i in range(len(mon_id_sorted)):
        ids.append("mon" + str(i+1) + " == \"" + mon_id_sorted[i][1] + "\"")

    my_nb_conditional_nograd = "    } else if (" + " and ".join(ids) + ") {\n"

    my_nb_conditional_nograd += "        mbnrg_{0}_deg{1}::mbnrg_{0}_deg{1}_{2} pot(".format(system_name,degree,version)

    ids = []
    for i in range(len(mon_id_sorted)):
        ids.append("mon" + str(mon_id_sorted[i][0] + 1))

    my_nb_conditional_nograd += ", ".join(ids) + ");\n"
    
    ids = []
    for i in range(len(mon_id_sorted)):
        ids.append("xyz" + str(mon_id_sorted[i][0] + 1) + ".data()")
    my_nb_conditional_nograd += "        energies = pot.eval(" + ", ".join(ids) + ", nm);\n"

    # Write code that needs to be added in the ONEBODY_GRD section of the code
    #cppout.write("=====>> SECTION ONEBODY_GRD <<=====\n")
    #cppout.write("File: src/potential/1b/energy1b.cpp\n")

    ids = []
    for i in range(len(mon_id_sorted)):
        ids.append("mon" + str(i+1) + " == \"" + mon_id_sorted[i][1] + "\"")

    my_nb_conditional_grad = "    } else if (" + " and ".join(ids) + ") {\n"

    my_nb_conditional_grad += "        mbnrg_{0}_deg{1}::mbnrg_{0}_deg{1}_{2} pot(".format(system_name,degree,version)

    ids = []
    for i in range(len(mon_id_sorted)):
        ids.append("mon" + str(mon_id_sorted[i][0] + 1))

    my_nb_conditional_grad += ", ".join(ids) + ");\n"
    
    ids = []
    idgs = []
    for i in range(len(mon_id_sorted)):
        ids.append("xyz" + str(mon_id_sorted[i][0] + 1) + ".data()")
        idgs.append("grad" + str(mon_id_sorted[i][0] + 1) + ".data()")
    my_nb_conditional_grad += "        energies = pot.eval(" + ", ".join(ids) + ", " + ", ".join(idgs) + " nm);\n"

    # Write code that needs to be added in the INCLUDE1B section of the code
    #cppout.write("=====>> SECTION INCLUDE1B <<=====\n")
    #cppout.write("File: src/potential/1b/energy1b.h\n")

    my_potential_include = "#include \"potential/" + str(number_of_monomers) + "b/mbnrg_{0}_deg{1}_{2}.h\"\n"




    print("Copying files to MBX_files folder...")

    # Folder where everything related to MBX is gonna go
    sofdir = workdir + "/MBX_files"

    # create folder with the files for the software
    os.system("mkdir -p " + sofdir)

    # define names for files
    headerf = "poly_{}b_{}_deg{}_{}.h".format(number_of_monomers,system_name,degree,version)
    cppgrad = "poly_{}b_{}_deg{}_grad_{}.cpp".format(number_of_monomers,system_name,degree,version)
    cppnograd = "poly_{}b_{}_deg{}_nograd_{}.cpp".format(number_of_monomers,system_name,degree,version)
    holderh = "mbnrg_{}b_{}_deg{}_{}.h".format(number_of_monomers,system_name,degree,version)
    holdercpp = "mbnrg_{}b_{}_deg{}_{}.cpp".format(number_of_monomers,system_name,degree,version)

    # Move them
    os.system("cp " + headerf + " " + cppgrad + " " + cppnograd + " " + holderh + " " + holdercpp + " " + sofdir)





    # Write the part of the code that needs to be put in dispersion
    #cppout.write("=====>> SECTION DISPERSION <<=====\n")
    #cppout.write("File: src/potential/dispersion2b.cpp\n")

    # atom_type_X will be like [0,0,0,1,1] for A3B2
    # atom_lable_X will be ['A','B'] for A3B2

    # Get the atomlist
    atom_types = []
    for fragment in monomers:
        atom_types.append(get_atom_types(fragment))

    # Get the list with actual types
    atom_types_number = []
    atom_types_letter = []
    for fragment in atom_types:
        atom_types_number.append([])
        atom_types_letter.append([])
        count = 0
        for i in range(1,len(fragment),2):
            for j in range(fragment[i]):
                atom_types_number[-1].append(count)
                atom_types_letter[-1].append(fragment[i-1])
            count += 1
     
    if number_of_monomers == 1:
        my_mon = [list(mon_id_sorted[0]), list(mon_id_sorted[0])]
        my_letter_types = [atom_types_letter[0],atom_types_letter[0]]
        my_number_types = [atom_types_number[0],atom_types_number[0]]
    elif number_of_monomers == 2:
        my_mon = [mon_id_sorted[0], mon_id_sorted[1]]
        my_letter_types = [atom_types_letter[0],atom_types_letter[1]]
        my_number_types = [atom_types_number[0],atom_types_number[1]]
    
    # First long range dispersion
    if number_of_monomers == 1:
        my_c6_lr_text = "    }} else if (mon_id == \"{}\") {{\n".format(mon_ids[0])
        my_c6_lr_text += "        for (size_t nv = 0; nv < n_mon; nv++) {{ \n"
        for j in range(len(atom_types_letter[0])):
            c6index = max(atom_types_number[0])*atom_types_number[0][j] + atom_types_number[0][j] 
            my_c6_long_range = C6[c6index]
            my_c6_lr_text += "            c6_lr[nv * natoms + fst_ind] = {}; // {}\n".format(math.sqrt(my_c6_long_range), atom_types_letter[0][j])


    if number_of_monomers < 3:
        if my_mon[0][0] == my_mon[1][0]:
            my_mon[1][0] += 1
        # Dispersion
        ids = []
        for i in range(len(my_mon)):
            ids.append("m" + str(i+1) + " == \"" + my_mon[i][1] + "\"")
        my_dispersion_text = "    } else if (" + " and ".join(ids) + ") {\n"
    
        for i in range(len(my_number_types)):
            for j in range(len(my_number_types[my_mon[i][0]])):
                my_dispersion_text += "        types{}.push_back({});\n".format(my_mon[i][0] + 1, my_number_types[my_mon[i][0]][j])
            my_dispersion_text += "\n"

        number_of_types_2 = 0
        if my_mon[0][0] > my_mon[1][0]:
            number_of_types_2 = max(my_number_types[my_mon[0][0]])
        else:
            number_of_types_2 = max(my_number_types[my_mon[1][0]])

        my_dispersion_text += "        nt2 = {};\n\n".format(number_of_types_2 + 1)

        my_dispersion_text += "        // Fill in (in order) the C6 and d6 coefficients\n"
    
        c6_units = "// kcal/mol * A^(-6) "
        d6_units = "// A^(-1)"
    
        c6_text = []
        d6_text = []
    
        for i in range(max(my_number_types[my_mon[0][0]]) + 1):
            for j in range(max(my_number_types[my_mon[1][0]]) + 1):
                c6index = max(my_number_types[my_mon[1][0]])*i + j
                let1 = my_letter_types[my_mon[0][0]][my_number_types[my_mon[0][0]].index(i)]
                let2 = my_letter_types[my_mon[1][0]][my_number_types[my_mon[1][0]].index(j)]
                c6_text.append("        C6.push_back(" + str(C6[c6index]) + ");  " + c6_units + " " + let1 + "--" + let2 + "\n")
                d6_text.append("        d6.push_back(" + str(d6[c6index]) + ");  " + d6_units + " " + let1 + "--" + let2 + "\n")

        for i in range(len(c6_text)):
            my_dispersion_text += c6_text[i]

        for i in range(len(d6_text)):
            my_dispersion_text += d6_text[i]

        # Repulsion
        ids = []
        for i in range(len(my_mon)):
            ids.append("m" + str(i+1) + " == \"" + my_mon[i][1] + "\"")
        my_buckingham_text = "    } else if (" + " and ".join(ids) + ") {\n"

        for i in range(len(my_number_types)):
            for j in range(len(my_number_types[my_mon[i][0]])):
                my_buckingham_text += "        types{}.push_back({});\n".format(my_mon[i][0] + 1, my_number_types[my_mon[i][0]][j])
            my_buckingham_text += "\n"

        number_of_types_2 = 0
        if my_mon[0][0] > my_mon[1][0]:
            number_of_types_2 = max(my_number_types[my_mon[0][0]])
        else:
            number_of_types_2 = max(my_number_types[my_mon[1][0]])

        my_buckingham_text += "        nt2 = {};\n\n".format(number_of_types_2 + 1)


        my_buckingham_text += "        // Fill in (in order) the C6 and d6 coefficients\n"

        A_units = "// kcal/mol"
        b_units = "// A^(-1)"

        A_text = []
        b_text = []

        shift = 0;
        for i in range(max(my_number_types[my_mon[0][0]]) + 1):
            for j in range(max(my_number_types[my_mon[1][0]]) + 1):
                c6index = max(my_number_types[my_mon[1][0]])*my_number_types[my_mon[0][0]][i] + my_number_types[my_mon[1][0]][j]
                let1 = my_letter_types[my_mon[0][0]][my_number_types[my_mon[0][0]].index(i)]
                let2 = my_letter_types[my_mon[1][0]][my_number_types[my_mon[1][0]].index(j)]
                A_text.append("        a.push_back(" + str(A_buck[c6index]) + ");  " + A_units + " " + let1 + "--" + let2 + "\n")
                b_text.append("        b.push_back(" + str(b_buck[c6index]) + ");  " + b_units + " " + let1 + "--" + let2 + "\n")

        for i in range(len(A_text)):
            my_buckingham_text += A_text[i]

        for i in range(len(b_text)):
            my_buckingham_text += b_text[i]


    # Write the C++ code in a file
    fcpp = open(sofdir + "/MBX_cpp_code.txt", 'w')

    # Print constructor information
    fcpp.write("// SECTION CONSTRUCTOR\n")
    fcpp.write("// " + holdercpp + "\n")
    fcpp.write(my_constructor_text)

    fcpp.write("\n\n\n// SECTION SITES\n")
    fcpp.write("src/bblock/sys_tools.h\n")
    fcpp.write(my_monomer_sites_text)

    fcpp.write("\n\n\n// SECTION CHARGES\n")
    fcpp.write("src/bblock/sys_tools.h\n")
    fcpp.write(my_monomer_charges_text)

    fcpp.write("\n\n\n// SECTION POLFACS\n")
    fcpp.write("src/bblock/sys_tools.h\n")
    fcpp.write(my_monomer_polfacs_text)

    fcpp.write("\n\n\n// SECTION POLS\n")
    fcpp.write("src/bblock/sys_tools.h\n")
    fcpp.write(my_monomer_pols_text)

    fcpp.write("\n\n\n// SECTION EXCLUDED\n")
    fcpp.write("src/bblock/sys_tools.h\n")
    fcpp.write(my_monomer_excluded_text)

    fcpp.write("\n\n\n SECTION INCLUDE{}B\n".format(number_of_monomers))
    fcpp.write("src/potential/{0}b/energy{0}b.h\n".format(number_of_monomers))
    fcpp.write(my_potential_include)

    fcpp.write("\n\n\n SECTION {}B_NO_GRADIENT\n".format(number_of_monomers))
    fcpp.write("src/potential/{0}b/energy{0}b.cpp\n".format(number_of_monomers))
    fcpp.write(my_nb_conditional_nograd)

    fcpp.write("\n\n\n SECTION {}B_GRADIENT\n".format(number_of_monomers))
    fcpp.write("src/potential/{0}b/energy{0}b.cpp\n".format(number_of_monomers))
    fcpp.write(my_nb_conditional_grad)

    fcpp.write("\n\n\n SECTION C6_LONG_RANGE\n")
    fcpp.write("src/bblock/sys_tools.h\n")
    fcpp.write(my_c6_lr_text)

    fcpp.write("\n\n\n SECTION DISPERSION\n")
    fcpp.write("src/bblock/sys_tools.h\n")
    fcpp.write(my_dispersion_text)

    fcpp.write("\n\n\n SECTION SITES\n")
    fcpp.write("src/bblock/sys_tools.h\n")
    fcpp.write(my_buckingham_text)

    fcpp.close()






