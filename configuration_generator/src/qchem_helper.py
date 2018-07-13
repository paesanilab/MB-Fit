import os
import sys
import numpy as np

#Gets the number of atoms by checking the first non-empty line
def get_natoms(xyzfile):
    with open(xyzfile, 'r') as f:
        for line in f:
            if line.strip():        
                natoms = line
                return int(natoms)

def get_coordinates(xyzfile):
    coordinates = """"""    
    with open(xyzfile, 'r') as f:
        coords_list = f.readlines()[2:]
    for coord in coords_list:
        coordinates += coord
    return coordinates
        
#Creates and writes an input file in the format for QChem
def write_qchem_input(file_name, charge, multiplicity, xyzfile, jobtype, method, basis, ecp):
    if ".inp" in file_name or ".in" in file_name:
        f = open(file_name, "w+")
    else:
        f = open(file_name+".inp", "w+")

    #$molecule section        
    f.write("$molecule\n")
    f.write("{} {}\n".format(charge, multiplicity))
    f.write(get_coordinates(xyzfile))
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
def optimize(infile_name, outfile_name, charge, multiplicity, xyzfile, method, basis, ecp):
    #Write the inputfile and call QChem    
    write_qchem_input(infile_name, charge, multiplicity, xyzfile, 'opt', method, basis, ecp)
    os.system("qchem %s %s" % (infile_name, outfile_name))
    
    found = False
    #found is set to true when the keyword 'Final energy is' is found. If found is true, it then looks for the keyword 'ATOM' to look for the optimized geometry
    parse = True
    molecule = []
    with open(outfile_name, "r") as outfile:
        for line in outfile:
            if parse:
                if found:
                    if "ATOM" in line:
                        for i in range(3):
                            molecule.append(next(outfile))
                        parse = False
                elif "Final energy is" in line:
                    e = line.split()[3]                
                    found = True
    return e, molecule
                
def frequencies(infile_name, outfile_name, charge, multiplicity, opt_geo_file, method, basis, ecp):
    write_qchem_input(infile_name, charge, multiplicity, opt_geo_file, 'freq', method, basis, ecp)
    os.system("qchem %s %s" % (infile_name, outfile_name))
    found = False
    parse = True
    num_atoms = get_natoms(opt_geo_file)
    modes_of_each_atom = []
    #modes_of_each_atom is a list that contains the mode of each atom line by line; it is not in the desired output format and is an intermediate step into getting the desired format
    normal_modes = []
    with open(outfile_name, "r") as outfile:
        for line in outfile:
            if parse:
                if found:
                    if "Mode:" in line:
                        num_modes = len(line.split())-1
                        frequencies = list(map(lambda x: float(x), next(outfile).split()[1:]))
                        next(outfile)
                        red_masses = list(map(lambda x: float(x), next(outfile).split()[2:]))
                        next(outfile)
                        next(outfile)
                        next(outfile)
                        next(outfile)
                        for atom in range(num_atoms):
                            modes_of_atom = list(filter(lambda x: x.strip(), next(outfile).split("   ")))[1:]
                            modes_of_atom = list(map(lambda x: x.split(), modes_of_atom))
                            for x in modes_of_atom:
                                for y in x:
                                    y = float(y)
                            modes_of_each_atom.append(modes_of_atom)  
                        for mode in range(num_modes):
                            array = np.zeros((num_atoms, 3))
                            for atom in range(num_atoms):
                                array[atom]=modes_of_each_atom[atom][mode]
                            normal_modes.append(array)
                        return normal_modes, frequencies, red_masses, num_atoms
                                       
                elif "INFRARED INTENSITIES (KM/MOL)" in line:
                    found = True
                
