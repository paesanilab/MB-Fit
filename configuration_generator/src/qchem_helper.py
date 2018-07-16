# qchem_helper
# @author April

import subprocess
import numpy as np
from molecule import Molecule
        
#Creates and writes an input file in the format for QChem
def write_qchem_input(file_name, charge, multiplicity, molecule, jobtype, method, basis, ecp):

    with open(file_name, 'w') as f:
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
        f.write("ecp " + ecp + "\n")
        f.write("$end\n")
        f.close()

#Uses the method defined above to create a QChem input file, then runs QChem, parses the file to look for the optimized geometry using keywords, and returns the energy as well as a list containing the lines of the optimized geometry in the output file
def optimize(molecule, filenames, config):
    #Write the inputfile and call QChem    
    write_qchem_input(filenames['qchem_opt_input'], config['molecule']['charges'], config['molecule']['spins'], molecule, 'opt', config['config_generator']['method'], config['config_generator']['basis'], config['config_generator']['ecp'])
    print("Optimizing geometry...")    
    subprocess.run("qchem %s %s" % (filenames['qchem_opt_input'], filenames['qchem_opt_output']), shell = True)
    
    found = False
    #found is set to true when the keyword 'Final energy is' is found. If found is true, it then looks for the keyword 'ATOM' to look for the optimized geometry
    qchem_mol = []
    with open(filenames['qchem_opt_output'], "r") as outfile:
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
    write_qchem_input(filenames['qchem_freq_input'], config['molecule']['charges'], config['molecule']['spins'], optimized_molecule, 'freq', config['config_generator']['method'], config['config_generator']['basis'], config['config_generator']['ecp'])
    print("Determining normal modes and running frequency analysis...")    
    subprocess.run("qchem %s %s" % (filenames['qchem_opt_input'], filenames['qchem_opt_output']), shell = True)
    found = False
    modes_of_each_atom = []
    #modes_of_each_atom is a list that contains the mode of each atom line by line; it is not in the desired output format and is an intermediate step into getting the desired format
    normal_modes = []
    with open(filenames['qchem_freq_output'], "r") as outfile:
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
                
