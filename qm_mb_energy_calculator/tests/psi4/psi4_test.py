try: #checks if psi4 present
    import psi4
    have_psi4 = True
except:
    have_psi4 = False

config_file = open("config.txt","r")

for line in config_file:

    line_parts = line.split()

    if "dev_output = " in line:
        dev_output = bool(line_parts[2])

    elif "nthreads = " in line:
        nthreads = int(line_parts[2])
    
    elif "initial_molecule_file_name = " in line:
        initial_molecule_file_name = line_parts[2]

    elif "method = " in line:
        method = line_parts[2]

    elif "basis = " in line:
        basis = line_parts[2]

    elif "psi4_output_file_name = " in line:
        psi4_output_file_name = line_parts[2]

    elif "normal_modes_output_file_name" in line:
        normal_modes_output_file = line_parts[2]

    elif "molecule_output_file_name" in line:
        molecule_output_file_name = line_parts[2]

    elif "charge = " in line:
        molecule_charge = line_parts[2]

    elif "multiplicity = " in line:
        molecule_multiplicity = line_parts[2]

    elif "psi4_memory_setting = " in line:
        psi4_memory_setting = line_parts[2]

config_file.close()

model = "{}/{}".format(method,basis)

#get molecule string
str_mol_rep = ""

initial_molecule_file = open(initial_molecule_file_name,"r")
for line in initial_molecule_file:
    str_mol_rep += line
initial_molecule_file.close()
str_mol_rep += "{} {}\n".format(molecule_charge, molecule_multiplicity)


if have_psi4:
    #Set psi4 settings
    psi4.core.set_output_file(psi4_output_file_name, False)
    psi4.set_memory(psi4_memory_setting)

    #Open Files
    out_file_dat = open(normal_modes_output_file,"w")
    in_file = open(psi4_output_file_name,"r")
    out_file_xyz = open(molecule_output_file_name,"w")

    #Normal Mode Variables Declaration and Instantiations
    correct_file_section = False
    count = 0
    atom_count = 0
    norm_modes = 0

    #START Regular Code :P
    mol = psi4.geometry(str_mol_rep)
    psi4.set_num_threads(nthreads)
    energy = psi4.energy(model) 
    opt = psi4.optimize(model, molecule = mol)
    e, wf = psi4.frequencies(model, molecule = mol, return_wfn=True)
    freq_array = wf.frequencies().to_array()

    #Development Output
    if dev_output:
        print(str_mol_rep)
        print("E({}) = {}".format(model, energy))
        print("")
        print("Frequencies: ")
        print(freq_array)
        print("")
        print("Normal Modes: ")

    for line in in_file: #traverses psi4.out file
        count += 1
        if ("Frequencies in cm^-1; force constants in au.") in line: #checks if on correct line of psi4.out file
            correct_file_section = True
            count = 0
        if correct_file_section and count > 4 and norm_modes < wf.frequencies().to_array().size:
            if atom_count < mol.natom(): #checks if in boundaries of one molecule's atoms, prints to out file and .bat file
                aspects = line.split("  ")

                if dev_output: #Dev Output
                    print(" " + aspects[1].rstrip() + " "  + aspects[2] + " "  + aspects[3] + " " + aspects[4].rstrip())

                if atom_count == 0: #checks if on first atom of current normal mode molecule
                    out_file_dat.write("Normal Mode:  {}\n".format(norm_modes + 1))
                    out_file_dat.write("{}\n".format(freq_array[norm_modes]))
                    out_file_dat.write("red_mass =          0\n")
                out_file_dat.write(" " * 9  + aspects[2] + " " * 9  + aspects[3] + " " * 9 + aspects[4].rstrip() + "\n")
                atom_count += 1
            else: #only if after one molecule finished, resets values and increments norm_modes
                atom_count = 0
                count = 1
                norm_modes += 1

                if dev_output: #Dev Output
                    print("")

    if dev_output: #Dev Output
        print("")
        print("String Representation of Molecule before optimization: ") #prints starting string representation of molecule
        print(str_mol_rep)
        print("")
        print("String Representation of Molecule after optimization: ") #prints final string representation of molecule
        print(mol.save_string_xyz())

    out_file_xyz.write("{}\n".format(mol.natom()))
    out_file_xyz.write("{}\n".format(e))
    out_file_xyz.write(mol.save_string_xyz())
    out_file_dat.close()
    out_file_xyz.close()
    in_file.close()
elif dev_output: #Dev Output
    print("psi4 not available")
