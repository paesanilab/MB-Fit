import os, sys
import json
import math
import filecmp
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

def generate_software_files(settings_path, config_file, mon_ids, do_ttmnrg, mbnrg_fits_path, degree, MBX_HOME, version,
                            virtual_sites):
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
        if any([virtual_site in monomers[i] for virtual_site in virtual_sites]):
            use_lonepairs[i] = 1
        if use_mbpol[i] != 0:
            number_of_sites[i] += 1
    
    # Obtain the lists with the excluded pairs
    excluded_pairs_12 = config.getlist("fitting", "excluded_pairs_12", int)
    excluded_pairs_13 = config.getlist("fitting", "excluded_pairs_13", int)
    excluded_pairs_14 = config.getlist("fitting", "excluded_pairs_14", int)
    vsites = config.getlist("fitting", "virtual_site_labels", str)
    
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
        A_buck = config.getlist("fitting", "a", float)
        b_buck = config.getlist("fitting", "d6", float)
        d6 = b_buck
    except:
        A_buck = [0.0] * len(C6)
        b_buck = [0.0] * len(C6)
        d6 = [0.0] * len(C6)
        print("***WARNING*** Seems like either A or d6 are not defined in the config file.These values will be filled as 0, but at this point of the process you should have the TTM-nrg fit perfermed.")

    ############################################################################
    ## Polynomial ##############################################################
    ############################################################################
    
    mon_id_sorted = sorted(list(enumerate(mon_ids)), key=lambda x: x[1])

    # Folder where everything related to MBX is gonna go
    workdir = os.getcwd()
    sofdir = workdir + "/MBX_files"

    # create folder with the files for the software
    os.system("mkdir -p " + sofdir)

    # Obtain polynomial coefficients and non_linear parameters
    
    if mbnrg_fits_path is not None: 
        print("Getting polynomial fitted parameters...")
        # Define a couple things
        mbnrg_best_fit = mbnrg_fits_path + "/best_fit"
        cdl_file = "fit-" + str(number_of_monomers) + "b.cdl"

        # Obtain polynomial coefficients and non_linear parameters
        constants = []
        polycoef = []
        npoly = -1
        
        with open(mbnrg_best_fit + "/" + cdl_file,'r') as cdl:
            line = cdl.readline()
            while True:
                # Skip name tag
                if line.strip().startswith(":"):
                    if not "name" in line:
                        constants.append(line.replace(":", "  m_"))
                # Find number of polynomial linear coefficients
                elif line.startswith("  poly"):
                    npoly = int(line.strip().split()[2].replace(";",""))
                # Store polynomial coefficients
                elif line.startswith("poly"):
                    for i in range(npoly):
                        polycoef.append("            " + cdl.readline().replace(";","};"))
                line = cdl.readline()
                if line == "":
                    break


        my_constructor_text = ""

        my_constructor_text += "    if ("
        ids = []
        for i in range(len(mon_id_sorted)):
            ids.append("mon" + str(i+1) + " == \"" + mon_ids[i] + "\"")

        my_constructor_text += " and ".join(ids) + ") {\n"

        my_constructor_text += "        coefficients = std::vector<double> {\n"

        ## Open file that contains the code to add
        #cppout = open("software_code.txt",'w')

        ## Write code that needs to be added in the constructor
        #cppout.write("=====>> SECTION CONSTRUCTOR <<=====\n")
        #cppout.write("File: newly generated x1b...degN_v1x.cpp\n")
        #a = '''
        for c in polycoef:
            my_constructor_text += c

        c += "\n\n"
        
        for p in constants:
            my_constructor_text += p

        my_constructor_text += "\n    } // end if " + " and ".join(ids) + "\n" 

    ############################################################################
    ## Monomer Properties ##### Only if one-body ###############################
    ############################################################################
    i = 100
    if number_of_monomers == 1 and mbnrg_fits_path is not None:
        print("Getting monomer properties...")
    # Write code that needs to be added in the SITES section of the code
        my_monomer_sites_text = ""
        my_monomer_sites_text += """
        }} else if (mon[i] == "{}") {{
            sites.push_back({});
            nat.push_back({});
""".format(mon_ids[0],number_of_sites[0],number_of_atoms[0])

    
    # Write code that needs to be added in the CHARGES section of the code
        my_monomer_charges_text = ""
        my_monomer_charges_text += """
        }} else if (mon_id == "{}") {{
            for (size_t nv = 0; nv < n_mon; nv++) {{
""".format(mon_ids[0])

        for j in range(len(charges[0])):
            my_monomer_charges_text += "                charges[fst_ind + nv*nsites + {}] = {} * CHARGECON;\n".format(j,charges[0][j])
        my_monomer_charges_text += "            }\n"


    # Write code that needs to be added in the POLFACS section of the code
        my_monomer_polfacs_text = ""
        my_monomer_polfacs_text += """
        }} else if (mon_id == "{}") {{
            for (size_t nv = 0; nv < n_mon; nv++) {{
""".format(mon_ids[0])

        for j in range(len(polarizability_factors[0])):
            my_monomer_polfacs_text += "                polfac[fst_ind + nv*nsites + {}] = {};\n".format(j,polarizability_factors[0][j])
        my_monomer_polfacs_text += "            }\n"


    # Write code that needs to be added in the POLS section of the code
        my_monomer_pols_text = ""
        my_monomer_pols_text += """
        }} else if (mon_id == "{}") {{
            for (size_t nv = 0; nv < n_mon; nv++) {{
""".format(mon_ids[0])

        for j in range(len(polarizabilities[0])):
            my_monomer_pols_text += "                pol[fst_ind + nv*nsites + {}] = {};\n".format(j,polarizabilities[0][j])
        my_monomer_pols_text += "            }\n"


    # Write code that needs to be added in excluded pair section
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

    if mbnrg_fits_path is not None:
        # Write code that needs to be added in the ONEBODY_NOGRD section of the code
        print("Getting energy calls...")

        system_name = "_".join(monomers)
        system_name = system_name.replace("(", "_o_").replace(")", "_c_")

        ids = []
        for i in range(len(mon_id_sorted)):
            ids.append("mon" + str(i+1) + " == \"" + mon_id_sorted[i][1] + "\"")

        my_nb_conditional_nograd = "    } else if (" + " and ".join(ids) + ") {\n"

        my_nb_conditional_nograd += "        mbnrg_{0}_deg{1}::mbnrg_{0}_deg{1}_{2} pot(".format(system_name,degree,version)

        mbx_order_to_poly_order = {}
        for mbx_index, mon_id in enumerate(mon_id_sorted):
            poly_index = [i[0] for i in mon_id_sorted].index(mbx_index)
            mbx_order_to_poly_order[mbx_index] = poly_index

        ids = []
        for i in range(len(mon_id_sorted)):
            ids.append("mon" + str(mbx_order_to_poly_order[i] + 1))

        my_nb_conditional_nograd += ", ".join(ids) + ");\n"
        
        ids = []
        if number_of_monomers == 1:
            return_key = "energies = "
        else:
            return_key = "return "
        for i in range(len(mon_id_sorted)):
            ids.append("xyz" + str(mbx_order_to_poly_order[i] + 1) + ".data()")
        my_nb_conditional_nograd += "        " + return_key + "pot.eval(" + ", ".join(ids) + ", nm);\n"

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
            ids.append("mon" + str(mbx_order_to_poly_order[i] + 1))

        my_nb_conditional_grad += ", ".join(ids) + ");\n"
        
        ids = []
        idgs = []
        if number_of_monomers == 1:
            return_key = "energies = "
        else:
            return_key = "energy = "
        for i in range(len(mon_id_sorted)):
            ids.append("xyz" + str(mbx_order_to_poly_order[i] + 1) + ".data()")
            idgs.append("grad" + str(mbx_order_to_poly_order[i] + 1) + ".data()")
        my_nb_conditional_grad += "        " + return_key + " pot.eval(" + ", ".join(ids) + ", " + ", ".join(idgs) + ", nm);\n"

        # Write code that needs to be added in the INCLUDE1B section of the code
        #cppout.write("=====>> SECTION INCLUDE1B <<=====\n")
        #cppout.write("File: src/potential/1b/energy1b.h\n")

        my_potential_include = "#include \"potential/" + str(number_of_monomers) + "b/mbnrg_{}b_{}_deg{}_{}.h\"\n".format(number_of_monomers,system_name,degree,version)

        # define names for files
        headerf = "poly_{}b_{}_deg{}_{}.h".format(number_of_monomers,system_name,degree,version)
        cppgrad = "poly_{}b_{}_deg{}_grad_{}.cpp".format(number_of_monomers,system_name,degree,version)
        cppnograd = "poly_{}b_{}_deg{}_nograd_{}.cpp".format(number_of_monomers,system_name,degree,version)
        holderh = "mbnrg_{}b_{}_deg{}_{}.h".format(number_of_monomers,system_name,degree,version)
        holdercpp = "mbnrg_{}b_{}_deg{}_{}.cpp".format(number_of_monomers,system_name,degree,version)

    # Write the part of the code that needs to be put in dispersion

    if number_of_monomers < 3:
    # Get the atomlist

        atom_types = []
        for fragment in monomers:
            atom_types.append(get_atom_types(fragment))

        print("Atom types:", atom_types)

        # Get the list with actual types
        atom_types_number = []
        atom_types_letter = []
        for fragment in atom_types:
            atom_types_number.append([])
            atom_types_letter.append([])
            count = 0

            atom_number_dict = {}

            for i in range(1,len(fragment),2):
                if fragment[i-1] in vsites:
                    continue
                try:
                    c = atom_number_dict[fragment[i-1]]
                except KeyError:
                    atom_number_dict[fragment[i-1]] = count
                    c = count
                    count += 1
                for j in range(fragment[i]):
                    atom_types_number[-1].append(c)
                    atom_types_letter[-1].append(fragment[i-1])

        print("Atom Types Number:", atom_types_number)
        print("Atom Types Letter:", atom_types_letter)
         
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
            my_c6_lr_text += "        for (size_t nv = 0; nv < n_mon; nv++) { \n"
            for j in range(len(atom_types_letter[0])):
                # formula for finding index of diagonal c6 constants in the list of c6 constants.
                # example: for A1B1C1, AA is at index 0, BB is at index 3, CC is at index 5
                c6index = sum(range(max(atom_types_number[0]) - atom_types_number[0][j] + 2, max(atom_types_number[0]) + 2))
                print(c6index)
                my_c6_long_range = C6[c6index]
                my_c6_lr_text += "            c6_lr[nv * natoms + fst_ind] = {}; // {}\n".format(math.sqrt(my_c6_long_range), atom_types_letter[0][j])
            my_c6_lr_text += "        }\n"

        if my_mon[0][0] == my_mon[1][0]:
            my_mon[1][0] += 1
        # Dispersion
        ids = []
        for i in range(len(my_mon)):
            ids.append("mon_id" + str(i+1) + " == \"" + my_mon[i][1] + "\"")
        my_dispersion_text = "    } else if (" + " and ".join(ids) + ") {\n"
    
        for i in range(len(my_number_types)):
            for j in range(len(my_number_types[my_mon[i][0]])):
                my_dispersion_text += "        types{}.push_back({});\n".format(i + 1, my_number_types[my_mon[i][0]][j])
            my_dispersion_text += "\n"

        number_of_types_2 = max(my_number_types[my_mon[1][0]])

        my_dispersion_text += "        nt2 = {};\n\n".format(number_of_types_2 + 1)

        my_dispersion_text += "        // Fill in (in order) the C6 and d6 coefficients\n"
    
        c6_units = "// kcal/mol * A^(-6) "
        d6_units = "// A^(-1)"

        c6_order_dict = {} # "A-B" -> 3
        c6_index = 0
        for let1 in my_letter_types[0]:
            for let2 in my_letter_types[1]:
                l1, l2 = sorted([let1, let2])
                try:
                    c6_order_dict[l1 + "-" + l2]
                except KeyError:
                    c6_order_dict[l1 + "-" + l2] = c6_index;
                    c6_index += 1
    
        c6_text = []
        d6_text = []
        for i in range(max(my_number_types[my_mon[0][0]]) + 1):
            for j in range(max(my_number_types[my_mon[1][0]]) + 1):
                let1, let2 = sorted([my_letter_types[my_mon[0][0]][my_number_types[my_mon[0][0]].index(i)],
                                    my_letter_types[my_mon[1][0]][my_number_types[my_mon[1][0]].index(j)]])
                c6index = c6_order_dict[let1 + "-" + let2]
                c6_text.append("        C6.push_back(" + str(C6[c6index]) + ");  " + c6_units + " " + let1 + "--" + let2 + "\n")
                d6_text.append("        d6.push_back(" + str(d6[c6index]) + ");  " + d6_units + " " + let1 + "--" + let2 + "\n")

        for i in range(len(c6_text)):
            my_dispersion_text += c6_text[i]

        for i in range(len(d6_text)):
            my_dispersion_text += d6_text[i]

        # Repulsion
        ids = []
        for i in range(len(my_mon)):
            ids.append("mon_id" + str(i+1) + " == \"" + my_mon[i][1] + "\"")
        my_buckingham_text = "    } else if (" + " and ".join(ids) + ") {\n"

        for i in range(len(my_number_types)):
            for j in range(len(my_number_types[my_mon[i][0]])):
                my_buckingham_text += "        types{}.push_back({});\n".format(i + 1, my_number_types[my_mon[i][0]][j])
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

        for i in range(max(my_number_types[my_mon[0][0]]) + 1):
            for j in range(max(my_number_types[my_mon[1][0]]) + 1):
                let1, let2 = sorted([my_letter_types[my_mon[0][0]][my_number_types[my_mon[0][0]].index(i)],
                                    my_letter_types[my_mon[1][0]][my_number_types[my_mon[1][0]].index(j)]])
                c6index = c6_order_dict[let1 + "-" + let2]
                A_text.append("        a.push_back(" + str(A_buck[c6index]) + ");  " + A_units + " " + let1 + "--" + let2 + "\n")
                b_text.append("        b.push_back(" + str(b_buck[c6index]) + ");  " + b_units + " " + let1 + "--" + let2 + "\n")

        for i in range(len(A_text)):
            my_buckingham_text += A_text[i]

        for i in range(len(b_text)):
            my_buckingham_text += b_text[i]


    print("Writing files to MBX_files folder...")

    # Write the C++ code in a file
    fcpp = open(sofdir + "/MBX_cpp_code.txt", 'w')

    
    if mbnrg_fits_path is not None:
        fcpp.write("// SECTION CONSTRUCTOR\n")
        fcpp.write("// " + holdercpp + "\n")
        fcpp.write(my_constructor_text)

        fcpp.write("\n\n\n SECTION INCLUDE{}B\n".format(number_of_monomers))
        fcpp.write("src/potential/{0}b/energy{0}b.h\n".format(number_of_monomers))
        fcpp.write(my_potential_include)

        fcpp.write("\n\n\n SECTION {}B_NO_GRADIENT\n".format(number_of_monomers))
        fcpp.write("src/potential/{0}b/energy{0}b.cpp\n".format(number_of_monomers))
        fcpp.write(my_nb_conditional_nograd)

        fcpp.write("\n\n\n SECTION {}B_GRADIENT\n".format(number_of_monomers))
        fcpp.write("src/potential/{0}b/energy{0}b.cpp\n".format(number_of_monomers))
        fcpp.write(my_nb_conditional_grad)

    if number_of_monomers == 1:
        fcpp.write("\n\n\n// SECTION SITES\n")
        fcpp.write("src/bblock/sys_tools.cpp\n")
        fcpp.write(my_monomer_sites_text)

        fcpp.write("\n\n\n// SECTION CHARGES\n")
        fcpp.write("src/bblock/sys_tools.cpp\n")
        fcpp.write(my_monomer_charges_text)

        fcpp.write("\n\n\n// SECTION POLFACS\n")
        fcpp.write("src/bblock/sys_tools.cpp\n")
        fcpp.write(my_monomer_polfacs_text)

        fcpp.write("\n\n\n// SECTION POLS\n")
        fcpp.write("src/bblock/sys_tools.cpp\n")
        fcpp.write(my_monomer_pols_text)

        fcpp.write("\n\n\n// SECTION EXCLUDED\n")
        fcpp.write("src/bblock/sys_tools.cpp\n")
        fcpp.write(my_monomer_excluded_text)

        fcpp.write("\n\n\n SECTION C6_LONG_RANGE\n")
        fcpp.write("src/bblock/sys_tools.cpp\n")
        fcpp.write(my_c6_lr_text)

    if do_ttmnrg:
        fcpp.write("\n\n\n SECTION DISPERSION\n")
        fcpp.write("src/potential/dispersion/disptools.cpp\n")
        fcpp.write(my_dispersion_text)

        fcpp.write("\n\n\n SECTION BUCKINGHAM\n")
        fcpp.write("src/potential/buckingham/bucktools.cpp\n")
        fcpp.write(my_buckingham_text)

    fcpp.close()

    if number_of_monomers == 1:
        this_mon = mon_ids[0]
        these_monomers = [mon_ids[0],mon_ids[0]]
    else:
        these_monomers = [mon_ids[0],mon_ids[1]]

    # Now we add the files to MBX if indicated
    if MBX_HOME is not None:
        # Start with monomer properties
        if number_of_monomers == 1:
            with open(MBX_HOME + "/src/bblock/sys_tools.cpp", 'r') as systools:
                lines = systools.readlines()
            # Sites
            i = 0
            while i < len(lines):
                if "BEGIN SECTION SITES" in lines[i]:
                    mon_there = False
                    while not "END SECTION SITES" in lines[i]:
                        if "\"" + this_mon + "\"" in lines[i]:
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_monomer_sites_text)
                    else:
                        print("***WARNING*** Seems like there is already a monomer {} defined in SECTION SITES. Won't replace anything.".format(this_mon))
                    break
                else:
                    i += 1
            # Charges
            i = 0
            while i < len(lines):
                if "BEGIN SECTION CHARGES" in lines[i]:
                    mon_there = False
                    while not "END SECTION CHARGES" in lines[i]:
                        if "\"" + this_mon + "\"" in lines[i]:
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_monomer_charges_text)
                    else:
                        print("***WARNING*** Seems like there is already a monomer {} defined in SECTION CHARGES. Won't replace anything.".format(this_mon))
                    break
                else:
                    i += 1
            # Polarizability
            i = 0
            while i < len(lines):
                if "BEGIN SECTION POLS" in lines[i]:
                    mon_there = False
                    while not "END SECTION POLS" in lines[i]:
                        if "\"" + this_mon + "\"" in lines[i]:
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_monomer_pols_text)
                    else:
                        print("***WARNING*** Seems like there is already a monomer {} defined in SECTION POLS. Won't replace anything.".format(this_mon))
                    break
                else:
                    i += 1
            # Polarizability Factor
            i = 0
            while i < len(lines):
                if "BEGIN SECTION POLFACS" in lines[i]:
                    mon_there = False
                    while not "END SECTION POLFACS" in lines[i]:
                        if "\"" + this_mon + "\"" in lines[i]:
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_monomer_polfacs_text)
                    else:
                        print("***WARNING*** Seems like there is already a monomer {} defined in SECTION POLFACS. Won't replace anything.".format(this_mon))
                    break
                else:
                    i += 1
            # Long Range C6
            i = 0
            while i < len(lines):
                if "BEGIN SECTION C6_LONG_RANGE" in lines[i]:
                    mon_there = False
                    while not "END SECTION C6_LONG_RANGE" in lines[i]:
                        if "\"" + this_mon + "\"" in lines[i]:
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_c6_lr_text)
                    else:
                        print("***WARNING*** Seems like there is already a monomer {} defined in SECTION C6_LONG_RANGE. Won't replace anything.".format(this_mon))
                    break
                else:
                    i += 1
            # Excluded pairs
            i = 0
            while i < len(lines):
                if "BEGIN SECTION EXCLUDED" in lines[i]:
                    mon_there = False
                    while not "END SECTION EXCLUDED" in lines[i]:
                        if "\"" + this_mon + "\"" in lines[i]:
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_monomer_excluded_text)
                    else:
                        print("***WARNING*** Seems like there is already a monomer {} defined in SECTION EXCLUDED. Won't replace anything.".format(this_mon))
                    break
                else:
                    i += 1
            with open(MBX_HOME + "/src/bblock/sys_tools.cpp", 'w') as systools:
                for line in lines:
                    systools.write(line)

        # Dispersion
        if number_of_monomers < 3:
            with open(MBX_HOME + "/src/potential/dispersion/disptools.cpp", 'r') as disptools:
                lines = disptools.readlines()
            i = 0
            while i < len(lines):
                if "BEGIN SECTION DISPERSION" in lines[i]:
                    mon_there = False
                    while not "END SECTION DISPERSION" in lines[i]:
                        if "\"" + these_monomers[0] + "\"" in lines[i] and "\"" + these_monomers[1] + "\"" in lines[i]:
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_dispersion_text)
                    else:
                        print("***WARNING*** Seems like there is already a dimer {} -- {} defined in SECTION DISPERSION. Won't replace anything.".format(these_monomers[0],these_monomers[1]))
                    break
                else:
                    i += 1

            with open(MBX_HOME + "/src/potential/dispersion/disptools.cpp", 'w') as disptools:
                for line in lines:
                    disptools.write(line)

            # Buckingham
            with open(MBX_HOME + "/src/potential/buckingham/bucktools.cpp", 'r') as bucktools:
                lines = bucktools.readlines()
            i = 0
            while i < len(lines):
                if "BEGIN SECTION BUCKINGHAM" in lines[i]:
                    mon_there = False
                    while not "END SECTION BUCKINGHAM" in lines[i]:
                        if "\"" + these_monomers[0] + "\"" in lines[i] and "\"" + these_monomers[1] + "\"" in lines[i]:
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_buckingham_text)
                    else:
                        print("***WARNING*** Seems like there is already a dimer {} -- {} defined in SECTION BUCKINGHAM. Won't replace anything.".format(these_monomers[0],these_monomers[1]))
                    break
                else:
                    i += 1

            with open(MBX_HOME + "/src/potential/buckingham/bucktools.cpp", 'w') as bucktools:
                for line in lines:
                    bucktools.write(line)


        if not ttm_only:
            # Include section
            with open(MBX_HOME + "/src/potential/{0}b/energy{0}b.h".format(number_of_monomers), 'r') as mbnrg:
                lines = mbnrg.readlines()
            i = 0
            while i < len(lines):
                if "BEGIN SECTION INCLUDE{}B".format(number_of_monomers) in lines[i]:
                    mon_there = False
                    while not "END SECTION INCLUDE{}B".format(number_of_monomers) in lines[i]:
                        if holderh in lines[i]:
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_potential_include)
                    else:
                        print("***WARNING*** Seems like there is already a {} file included in SECTION INCLUDE{}B. Won't replace anything.".format(holderh, number_of_monomers))
                    break
                else:
                    i += 1

            with open(MBX_HOME + "/src/potential/{0}b/energy{0}b.h".format(number_of_monomers), 'w') as mbnrg:
                for line in lines:
                    mbnrg.write(line)

            # Energy with no gradients section
            with open(MBX_HOME + "/src/potential/{0}b/energy{0}b.cpp".format(number_of_monomers), 'r') as mbnrg:
                lines = mbnrg.readlines()
            i = 0
            while i < len(lines):
                if "BEGIN SECTION {}B_NO_GRADIENT".format(number_of_monomers) in lines[i]:
                    mon_there = False
                    while not "END SECTION {}B_NO_GRADIENT".format(number_of_monomers) in lines[i]:
                        if all("\"" + x + "\"" in lines[i] for x in mon_ids):
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_nb_conditional_nograd)
                    else:
                        nmer = " -- ".join(mon_ids)
                        print("***WARNING*** Seems like there is already a n-mer {} defined in SECTION {}B_NO_GRADIENT. Won't replace anything.".format(nmer, number_of_monomers))
                    break
                else:
                    i += 1

            # Energy with gradients section
            i = 0
            while i < len(lines):
                if "BEGIN SECTION {}B_GRADIENT".format(number_of_monomers) in lines[i]:
                    mon_there = False
                    while not "END SECTION {}B_GRADIENT".format(number_of_monomers) in lines[i]:
                        if all("\"" + x + "\"" in lines[i] for x in mon_ids):
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_nb_conditional_grad)
                    else:
                        nmer = " -- ".join(mon_ids)
                        print("***WARNING*** Seems like there is already a n-mer {} defined in SECTION {}B_GRADIENT. Won't replace anything.".format(nmer, number_of_monomers))
                    break
                else:
                    i += 1

            with open(MBX_HOME + "/src/potential/{0}b/energy{0}b.cpp".format(number_of_monomers), 'w') as mbnrg:
                for line in lines:
                    mbnrg.write(line)

            # Update constructor
            mbnrg_exists = os.path.exists(MBX_HOME + "/src/potential/{0}b/".format(number_of_monomers) + holdercpp)
            # If the holder exists, we need to add the constructor info there
            # If not, we will update the one in MBX_files
            if mbnrg_exists:
                # In this case, let's make sure that the polynomial holders and 
                # polynomials are actually the same
                # Compare polynomial with grads
                we_are_good = filecmp.cmp(MBX_HOME + "/src/potential/{0}b/".format(number_of_monomers) + cppgrad, sofdir + "/" + cppgrad)
                # Compare polynomial with no grads
                we_are_good = we_are_good and filecmp.cmp(MBX_HOME + "/src/potential/{0}b/".format(number_of_monomers) + cppnograd, sofdir + "/" + cppnograd)
                # Compare polynomial header
                we_are_good = we_are_good and filecmp.cmp(MBX_HOME + "/src/potential/{0}b/".format(number_of_monomers) + headerf, sofdir + "/" + headerf)
                # poly holder header
                we_are_good = we_are_good and filecmp.cmp(MBX_HOME + "/src/potential/{0}b/".format(number_of_monomers) + holderh, sofdir + "/" + holderh)
                if not we_are_good:
                    print("***WARNING*** Although the file {} already exists in {}/src/potential/{}b/, the polynomial files do not seem to match with the ones generated here. Parameters will be added in {}, but nothing will be copied to MBX for safety. Do that manually.".format(holdercpp, MBX_HOME, number_of_monomers, sofdir))
            else:
                we_are_good = True

            copy_holder = False
            if mbnrg_exists and we_are_good:
                mbnrg_file = MBX_HOME + "/src/potential/{0}b/".format(number_of_monomers) + holdercpp
            else:
                mbnrg_file = sofdir + "/" + holdercpp
                copy_holder = True
            with open(mbnrg_file, 'r') as mbnrg:
                lines = mbnrg.readlines()
            i = 0
            while i < len(lines):
                if "BEGIN SECTION CONSTRUCTOR" in lines[i]:
                    mon_there = False
                    while not "END SECTION CONSTRUCTOR" in lines[i]:
                        if all("\"" + x + "\"" in lines[i] for x in mon_ids):
                            mon_there = True
                            break
                        i += 1
                    if not mon_there:
                        lines.insert(i,my_constructor_text)
                    else:
                        nmer = " -- ".join(mon_ids)
                        print("***WARNING*** Seems like there is already a n-mer {} defined in SECTION CONSTRUCTOR. Won't replace anything.".format(nmer))
                    break
                else:
                    i += 1

            with open(mbnrg_file, 'w') as mbnrg:
                for line in lines:
                    mbnrg.write(line)

            if we_are_good:
                os.system("cd " + sofdir + " ; cp " + headerf + " " + cppgrad + " " + cppnograd + " " + holderh + " " + MBX_HOME + "/src/potential/{0}b/".format(number_of_monomers) + "; cd ../ ")

                # Now update CMake
                with open(MBX_HOME + "/src/potential/{0}b/".format(number_of_monomers) + "CMakeLists.txt", 'r') as cmakelists:
                    lines = cmakelists.readlines()
                i = 0
                while i < len(lines):
                    if "BEGIN SECTION CMAKELISTS" in lines[i]:
                        holder_there = False
                        polygrad_there = False
                        polynograd_there = False
                        while not "END SECTION CMAKELISTS" in lines[i]:
                            if cppgrad in lines[i]:
                                polygrad_there = True
                            if cppnograd in lines[i]:
                                polynograd_there = True
                            if holdercpp in lines[i]:
                                holder_there =  True
                            something_there = holder_there or polynograd_there or polygrad_there
                            if something_there:
                                break
                            i += 1
                        if not something_there:
                            lines.insert(i,"                 " + cppgrad + "\n")
                            lines.insert(i,"                 " + cppnograd + "\n")
                            lines.insert(i,"                 " + holdercpp + "\n")
                        else:
                            nmer = " or ".join([cppgrad, cppnograd, holdercpp])
                            print("***WARNING*** Seems like there are already some of these files ({}) in CMakeLists... Won't replace anything.".format(nmer))
                        break
                    else:
                        i += 1

                with open(MBX_HOME + "/src/potential/{0}b/".format(number_of_monomers) + "CMakeLists.txt", 'w') as cmakelists:
                    for line in lines:
                        cmakelists.write(line)


                if copy_holder:
                    os.system("cp " + sofdir + "/" + holdercpp + " " + MBX_HOME + "/src/potential/{0}b/".format(number_of_monomers))
















