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
    subprocess.run("qchem -nt %d %s %s" % (nt, fninp, fnout), shell = True)
    
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
    subprocess.run("qchem -nt %d %s %s" % (nt, fninp, fnout), shell = True)

    found = False
    modes_of_each_atom = []
    #modes_of_each_atom is a list that contains the mode of each atom line by line; it is not in the desired output format and is an intermediate step into getting the desired format
    normal_modes = []
    with open(fnout, "r") as outfile:
        for line in outfile:
            if found:
                if "Mode:" in line:
                    #Extracts the number of modes from the line with 'Mode: 1, 2, 3..'                        
                    num_modes = len(line.split())-1
                    #Extracts the frequencies from the line after
                    frequencies = list(map(lambda x: float(x), next(outfile).split()[1:]))
                    next(outfile)
                    #Extracts the reduced masses from the line 2 lines after
                    red_masses = list(map(lambda x: float(x), next(outfile).split()[2:]))
                    for i in range(4):
                        next(outfile)
                    #Gets the normal modes from the lines 4 lines after
                    for atom in range(optimized_molecule.num_atoms):
                        #Each line is split into a list modes_of_atom, containing the x y z coordinate of each normal mode of the element in order.
                        #strip() removes all the strings containing only spaces in the .split() list. The [1:] removes the element name and returns only the normal modes
                        modes_of_atom = list(filter(lambda x: x.strip(), next(outfile).split("   ")))[1:]
                        #Each normal mode, originally in a string, is split into a list [x, y, z]
                        modes_of_atom = list(map(lambda x: x.split(), modes_of_atom))
                        #Converts the normal mode coordinates from string to float
                        for x in modes_of_atom:
                            for y in x:
                                y = float(y)
                        #This list of the normal modes of 1 atom, modes_of_atom, is appended to a larger list modes_of_each_atom
                        modes_of_each_atom.append(modes_of_atom) 
                    #Creates an empty array and then fills it in to list out each normal mode 
                    for mode in range(num_modes):
                        array = np.zeros((optimized_molecule.num_atoms, 3))
                        for atom in range(optimized_molecule.num_atoms):
                            array[atom]=modes_of_each_atom[atom][mode]
                        normal_modes.append(array)
                    return normal_modes, frequencies, red_masses
                 
                #Uses "INFRARED INTENSITIES (KM/MOL)" as a keyword to find where the normal modes are                      
            elif "INFRARED INTENSITIES (KM/MOL)" in line:
                found = True
                
