# qchem_helper
# @author April

import subprocess
import numpy as np
from molecule import Molecule
        
#Creates and writes an input file in the format for QChem
def write_qchem_input(fn, config, molecule, jobtype):

    charge = config['molecule']['charges']        #AWG only works for monomer
    multiplicity = config['molecule']['spins']    #AWG -"-
    method = config['config_generator']['method']
    basis = config['config_generator']['basis']
    try:
        ecp = config['config_generator']['ecp']
    except:
        ecp = None

    with open(fn, 'w') as f:
        #$molecule section        
        f.write("$molecule\n")
        f.write("{} {}\n".format(charge, multiplicity))
        f.write(str(molecule))
        f.write("$end\n")
        f.write("\n")

        #$rem variables
        f.write("$rem\n")
        f.write("jobtype " + jobtype + "\n")
        f.write("method " + method + "\n")
        f.write("basis " + basis + "\n")
        if ecp:
            f.write("ecp " + ecp + "\n")
        f.write("$end\n")
        f.close()

# 1 Create a QChem input file
# 2 run QChem,
# 3 parse the file to look for the optimized geometry using keywords
# 4 return the energy and a list containing the lines of the optimized geometry
def optimize(molecule, filenames, config):

    #Write the inputfile and call QChem    
    fninp = filenames["log_name"] + "_opt.inp"
    fnout = filenames["log_name"] + "_opt.out"
    write_qchem_input(fninp, config, molecule, 'opt')

    print("Optimizing geometry...") 
    try:
        nt = int(config["qchem"]["num_threads"])
    except:
        nt = 1
    subprocess.run("qchem -nt %d %s %s" % (nt, fninp, fnout), 
                   stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                   shell=True, check=True)
    
    found = False
    #found is set to true when the keyword 'Final energy is' is found. If found is true, it then looks for the keyword 'ATOM' to look for the optimized geometry
    qchem_mol = []
    with open(fnout, "r") as outfile:
        for line in outfile:
            if found:
                if "ATOM" in line:
                    for i in range(molecule.num_atoms):
                        qchem_mol.append(next(outfile))
                        #Returns a list containing the lines with the optimized geometry coordinates
                    return float(e), qchem_mol
            elif "Final energy is" in line:
                e = line.split()[3]                
                found = True

def read_qchem_mol(qchem_mol, num_atoms, au_conversion = 1.0):
    molecule = """\n\n"""
    for line in qchem_mol:
        molecule += line[9:]
    return Molecule(molecule, au_conversion = 1.0)

def frequencies(optimized_molecule, filenames, config):

    fninp = filenames["log_name"] + "_freq.inp"
    fnout = filenames["log_name"] + "_freq.out"
    write_qchem_input(fninp, config, optimized_molecule, 'freq')

    print("Determining normal modes and running frequency analysis...")    
    try:
        nt = int(config["qchem"]["num_threads"])
    except:
        nt = 1
    subprocess.run("qchem -nt %d %s %s" % (nt, fninp, fnout), 
                   stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                   shell=True, check=True)

    with open(fnout, "r") as qchem_out:
        return read_normal_modes(qchem_out, optimized_molecule.num_atoms)


def read_normal_modes(file, num_atoms):
    frequencies = []
    red_masses = []
    normal_modes = []
    while(True):
        line = file.readline()

        if line == "":
            break

        if "Mode:" in line:
            num_modes = len(line.split()) - 1

            # read frequencies from file and cast each to a float
            frequencies += [float(freq) for freq in file.readline().split()[1:]]

            # skip the Force Constant line
            file.readline()

            # read red masses from file and cast each to float
            red_masses += [float(freq) for freq in file.readline().split()[2:]]

            # skip IR Active, IR Intens, Raman Active and labels line
            for i in range(4):
                file.readline()

            # initialize 3d list for normal modes
            modes = [[[0, 0, 0] for i in range(num_atoms)] for k in range(num_modes)]

            for atom in range(num_atoms):

                # read the coordinates for the next atom for the next num_modes modes
                modes_of_atom = [float(i) for i in file.readline().split()[1:]]

                # loop over each mode
                for index in range(num_modes):

                    # store coordinates in corresponding entries of modes list
                    modes[index][atom][0] = modes_of_atom[index * 3]
                    modes[index][atom][1] = modes_of_atom[index * 3 + 1]
                    modes[index][atom][2] = modes_of_atom[index * 3 + 2]

            normal_modes += modes

    return normal_modes, frequencies, red_masses
