import os, sys
import json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../../")
import configparser

# TODO Duplicated code. Maybe it can be imported from somwhere else?
def parse_array(string):
    num_open_brackets = 0
    elements = []
    element = ""
    for character in string:
        if num_open_brackets == 1 and character == ",":
            elements.append(element)
            element = ""
        elif character == "[":
            if num_open_brackets != 0:
                element += "["
            num_open_brackets += 1
        elif character == "]":
            num_open_brackets -= 1
            if num_open_brackets != 0:
                element += "]"
        else:
            element += character

    if num_open_brackets == 0:
        elements.append(element)
    else:
        print("Something went wrong while parsing a string into a list")
    return [parse_array(element) if "," in element else element for element in elements]


def write_cpp_software_properties(settings, config_file, fit_path, fit_cdl, mon_name1, mon_name2, poly_order, molecules):
    # NOTE mon_name is the actual monomer name (h2o, ch4, c2h6...)

    # read config
    config = configparser.ConfigParser()
    config.read(settings)
    config.read(config_file)

    vsites = ["X","Y","Z"]

    # Number of atoms in monomer
    # Store them in nat1 and nat2
    nat1, nat2 = (int(num_atoms) for num_atoms in config["molecule"]["fragments"].split(","))
#    nat = parse_array(config["molecule"]["fragments"])
#    nat1, nat2 = (int(num_atoms) for num_atoms in nat)

    # Symmetry of the molecules
    mol1 = molecules[0]
    mol2 = molecules[1]

    # Get types of atoms in each mol, needed for dispersion
    types_a = list(mol1)
    types_b = list(mol2)

    # Get list of types needed later
    # atom_type_X will be like [0,0,0,1,1] for A3B2
    # atom_lable_X will be ['A','B'] for A3B2
    atom_type_a = []
    atom_label_a = []
    t = 0

    for type_index in range(0, len(types_a), 2):
        if types_a[type_index] in vsites:
            continue

        for atom_index in range(1, int(types_a[type_index + 1]) + 1):
            atom_type_a.append(t)
        t += 1
        atom_label_a.append(types_a[type_index])
    
    atom_type_b = []
    atom_label_b = []
    t = 0
    
    for type_index in range(0, len(types_b), 2):
        if types_b[type_index] in vsites:
            continue

        for atom_index in range(1, int(types_b[type_index + 1]) + 1):
            atom_type_b.append(t)
        t += 1
        atom_label_b.append(types_b[type_index])
    
    # Monomer names. Ordering from lower to higher if necessary
    # co2 < h2o , h2o < rb
    sorted_mon = 0

    if mon_name1 > mon_name2:
        mon1 = mon_name2
        mon2 = mon_name1
        sorted_mon = 1
    else:
        mon1 = mon_name1
        mon2 = mon_name2
    
    # Get C6 from config
    c6_constants = parse_array(config["fitting"]["c6"])
    C6 = c6_constants[len(c6_constants) - 1]
    # And d6
    d6_constants = parse_array(config["fitting"]["d6"])
    d6 = d6_constants[len(d6_constants) - 1]

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
    if (mon1 == "''' + mon1 + '''" && mon2  == "''' + mon2 + '''") {
        coefficients = std::vector<double> { 
'''
    cppout.write(a)
    for c in polycoef:
        cppout.write(c)

    cppout.write("\n")

    for p in constants:
        cppout.write(p)

    a = '''
    } // if mon1 == "''' + mon1 + '''" && mon2  == "''' + mon2 + '''"




'''

    cppout.write(a)

    # Write code that needs to be added in the ONEBODY_NOGRD section of the code
    cppout.write("=====>> SECTION TWOBODY_NOGRD <<=====\n")
    cppout.write("File: src/potential/2b/energy2b.cpp\n")

    if sorted_mon == 1:
        a = '''
    } else if (m1 == "''' + mon1 + '''" && m2 == "''' + mon2 + '''") {
        x2b_''' + mol1 + "_" + mol2 + "_deg" + str(poly_order) + '''::x2b_''' + mol1 + "_" + mol2 + '''_v1x pot(m1,m2);
        return pot.eval(xyz2.data(), xyz1.data(), nm);
'''
    else:
        a = '''
    } else if (m1 == "''' + mon1 + '''" && m2 == "''' + mon2 + '''") {
        x2b_''' + mol1 + "_" + mol2 + "_deg" + str(poly_order) + '''::x2b_''' + mol1 + "_" + mol2 + '''_v1x pot(m1,m2);
        return pot.eval(xyz1.data(), xyz2.data(), nm);


'''

    cppout.write(a)

    # Write code that needs to be added in the ONEBODY_GRD section of the code
    cppout.write("=====>> SECTION TWOBODY_GRD <<=====\n")
    cppout.write("File: src/potential/2b/energy2b.cpp\n")

    if sorted_mon == 1:
        a = '''
    } else if (m1 == "''' + mon1 + '''" && m2 == "''' + mon2 + '''") {
        x2b_''' + mol1 + "_" + mol2 + "_deg" + str(poly_order) + '''::x2b_''' + mol1 + "_" + mol2 + '''_v1x pot(m1,m2);
        energy = pot.eval(xyz2.data(), xyz1.data(), grd2.data(), grd1.data(), nm);


'''
    else:
        a = '''
    } else if (m1 == "''' + mon1 + '''" && m2 == "''' + mon2 + '''") {
        x2b_''' + mol1 + "_" + mol2 + "_deg" + str(poly_order) + '''::x2b_''' + mol1 + "_" + mol2 + '''_v1x pot(m1,m2);
        energy = pot.eval(xyz1.data(), xyz2.data(), grd1.data(), grd2.data(), nm);


'''

    cppout.write(a)


    # Write code that needs to be added in the INCLUDE1B section of the code
    cppout.write("=====>> SECTION INCLUDE2B <<=====\n")
    cppout.write("File: src/potential/2b/energy2b.h\n")

    cppout.write("#include \"potential/2b/x2b_"  + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x.h\"\n\n\n\n")

    
    # Write the part of the code that needs to be put in dispersion
    cppout.write("=====>> SECTION DISPERSION <<=====\n")
    cppout.write("File: src/potential/dispersion2b.cpp\n")


    if sorted_mon == 0:
        a = '''
    } else if (m1 == "''' + mon1 + '''" and m2 == "''' + mon2 + '''") {
        // Define number of atoms in each mon
        nat1 = ''' + str(nat1) + ''';
        nat2 = ''' + str(nat2) + ''';

        // Define the type of atom in each mon
'''
    
        cppout.write(a)
    
        for i in atom_type_a:
            cppout.write("        types1.push_back(" + str(i) + ");\n")
    
        cppout.write("\n")
    
        for i in atom_type_b:
            cppout.write("        types2.push_back(" + str(i) + ");\n")
    
        cppout.write("\n")
    
        a = '''
        // Set the number of different types
        nt2 = ''' + str(len(atom_label_b)) + ''';

        // Fill in (in order) the C6 and d6 coefficients
'''
    
        cppout.write(a)
    
        c6_units = "// kcal/mol * A^(-6) "
        d6_units = "// A^(-1)"
    
        c6_text = []
        d6_text = []
    
        shift = 0;
        for i in range(len(atom_label_a)):
            if mol1 == mol2:
                shift -= i 
    
            for j in range(len(atom_label_b)):
                c6index = shift + len(atom_label_b)*i + j
                
                c6_text.append("        C6.push_back(" + C6[c6index] + ");  " + c6_units + " " + atom_label_a[i] + "--" + atom_label_b[j] + "\n")
                d6_text.append("        d6.push_back(" + d6[c6index] + ");  " + d6_units + " " + atom_label_a[i] + "--" + atom_label_b[j] + "\n")

    else:
        a = '''
    } else if (m1 == "''' + mon1 + '''" and m2 == "''' + mon2 + '''") {
        // Define number of atoms in each mon
        nat2 = ''' + str(nat1) + ''';
        nat1 = ''' + str(nat2) + ''';

        // Define the type of atom in each mon
'''

        cppout.write(a)

        for i in atom_type_a:
            cppout.write("        types2.push_back(" + str(i) + ");\n")

        cppout.write("\n")

        for i in atom_type_b:
            cppout.write("        types1.push_back(" + str(i) + ");\n")

        cppout.write("\n")

        a = '''
        // Set the number of different types
        nt2 = ''' + str(len(atom_label_a)) + ''';

        // Fill in (in order) the C6 and d6 coefficients
'''

        cppout.write(a)

        c6_units = "// kcal/mol * A^(-6) "
        d6_units = "// A^(-1)"

        c6_text = []
        d6_text = []

        shift = 0;
        for j in range(len(atom_label_b)):
            for i in range(len(atom_label_a)):
                c6index = shift + len(atom_label_b)*i + j

                c6_text.append("        C6.push_back(" + C6[c6index] + ");  " + c6_units + " " + atom_label_b[j] + "--" + atom_label_a[i] + "\n")
                d6_text.append("        d6.push_back(" + d6[c6index] + ");  " + d6_units + " " + atom_label_b[j] + "--" + atom_label_a[i] + "\n")

    for i in c6_text:
        cppout.write(i)

    cppout.write("\n")

    for i in d6_text:
        cppout.write(i)

    cppout.write("\n\n\n\n")

    cppout.close()

    # software folder name
    sofdir = fit_path + "/" + "software_files"

    os.system("mv software_code.txt " + sofdir)


def generate_software_files_2b(settings, symmetry, poly_path, poly_order, fit_path, config_file, fit_cdl, mon_name1, mon_name2):
    """
    Generates the parts of the C++ code needed to add the PEF to the energy software

    Args:
        settings - the file containing all relevent settings information
        symmetry - the A3B2.in type file
        poly_path   - directory where polynomial files are
        poly_order - the order of the polynomial in poly_path
        fit_path - directory to generate fit code in
        config_file    - monomer config file
        fit_cdl - the output cdl file of the fit to use
        mon_name - the human understandable name of the monomer (co2, h2o, no3...)

    Returns:
        None
    """

    # software folder name
    sofdir = fit_path + "/" + "software_files"

    # create folder with the files for the software
    os.system("mkdir -p " + sofdir)

    # get just the two A2B3 and C1D2 parts of input.in
    molecules = symmetry.split("_")
    mol1 = molecules[0]
    mol2 = molecules[1]

    # Create mon_names list and check if we reorder them:
    sorted_mon = 0

    if mon_name1 > mon_name2:
        mon1 = mon_name2
        mon2 = mon_name1
        sorted_mon = 1
    else:
        mon1 = mon_name1
        mon2 = mon_name2

    
    
    # mv x1b files needed
    os.system("mv x2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x.cpp " + sofdir)
    os.system("mv x2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x.h " + sofdir)

    # define names for files
    headerf = sofdir + "/poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x.h"
    cppgrad = sofdir + "/poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x.cpp"
    cppnograd = sofdir + "/poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1.cpp"

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
        newline = lines[line_num].replace("POLY_MODEL","POLY_2B_" + mol1 + "_" + mol2 + "_DEG" + str(poly_order)).replace("mb_system","x2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order)).replace("poly_model","poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x").replace("poly-model","poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x")

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
        newline = lines[line_num].replace("poly_2b_" + mol1 + "_" + mol2 + "_v1x.h","poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order)).replace("mb_system","x2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order)).replace("poly_model","poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x").replace("poly-model","poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x")
        
        line_num += 1
        grd.write(newline)

    # modify cpp nograd file
    with open(cppnograd + ".tmp", "r") as grd:
        lines = grd.readlines()

    grd = open(cppnograd, "w")
    line_num = 0
    while line_num < len(lines):
        # go line by line, and replace keywords
        newline = lines[line_num].replace("poly_2b_" + mol1 + "_" + mol2 + "_v1x.h","poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order)).replace("mb_system","x2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order)).replace("poly_model","poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x").replace("poly-model","poly_2b_" + mol1 + "_" + mol2 + "_deg" + str(poly_order) + "_v1x")

        line_num += 1
        grd.write(newline)

    grd.close()

    # remove tmp files
    os.system("rm " + cppgrad + ".tmp " + cppnograd + ".tmp " + headerf + ".tmp")

    # Call the function to write the individual pieces of code, defined above
    write_cpp_software_properties(settings, config_file, fit_path, fit_cdl, mon_name1, mon_name2, poly_order, molecules)

