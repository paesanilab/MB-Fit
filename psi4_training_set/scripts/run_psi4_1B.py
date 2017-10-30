#!/usr/bin/env python3

import os, sys
import psi4
import pybel as pb
import shutil
import argparse

from reader import readfile

#TODO
#Things to add later
#descriptive names for variables etc.
#reader for job control settings for psi4 like charge and multiplicity
#output storage and compression
#training set writing during/after calculation

# TODO (done)
# Find a way to parse the test input file
# Cannot use BOTH pybel and psi4
# The format for xyz input file:
# amount of atoms in a molecule
# comment(mandatory, but can be empty)
# n atoms, and their xyz coordinates

# Error flags, caused by invalid command-line options
memErr = methodErr = basisErr = fileErr = 0

'''
Subclass and function to override argparse
'''
class optionParse(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_lines):
        return arg_lines.split()[1:]

# Check for initial passed-in command-line args
parser = optionParse(fromfile_prefix_chars='@')


parser.add_argument('memory', type=str)
parser.add_argument('method', type=str) 
parser.add_argument('basis', type=str)
parser.add_argument('charge', type=int)
parser.add_argument('spin', type=int)
parser.add_argument('convergence', type=int)
parser.add_argument('cycles', type=int)
parser.add_argument('threshold', type=int)
parser.add_argument('input', type=str)

args = parser.parse_args()
print(args)

# If only one argument (meaning no arguments passed in),
# use default arguments
# If all arguments are passed in, check if arguments are valid(?)
# Else, do not continue

training_set_file = open('./training_set.xyz', 'w')

#read in an xyz file called input.xyz and start counting 
#configurations with i_config


molCount = 0

# This section is currently defunct. Will uncomment when issue is resolved
'''
for mol in pb.readfile("xyz","input.xyz"):
    i_config += 1

    #convert "mol" to string 
    my_mol = mol.write("xyz")

    #count the number of atoms per "mol" and create a range for iteration
    atom_total_range = range(len(mol.atoms))  

    #create the xyz_string to be used in psi4
    xyz_string = list()
    for atom in atom_total_range:
    
        xyz_string.append(my_mol.splitlines()[atom+2])
 
    xyz_string = "\n".join(xyz_string)
'''

    #TODO
    #output file necessary! we just need the energies of the configurations 
    #for now the output files will be placed somewhere but 
    #right now the script creates a directory called "calculations"
    #also creates a subdirectory for the output of each configuration 
    #with hardcoded 3 leading zeroes for the subdirectory name

# Temporary solution to pybel conflict. See reader.py
try:
    molecules = readfile('input.xyz')
except:
    print("Invalid file. Exiting.")
    training_set_file.close()

# Perform a calculation for each molecule
for mol in molecules:
    molCount += 1
    molString = mol.toString()

    if not os.path.exists('./calculations/'+str(molCount).zfill(4)):
        os.makedirs('./calculations/'+str(molCount).zfill(4))

    psi4.core.set_output_file('./calculations/'+str(molCount).zfill(4)+
        '/output.dat', False)
    
    #! Sample HF/3-21g one-body Computation
    # There are several ways we can do this without hard-coding.
    # 1. Use command-line input (user-defined), and set a default value if no
    #    input given
    # 2. Figure out how much memory do we actually need for a calculation,
    #    based on the molecule complexity
    try:
      psi4.set_memory(args.memory)
    except:
      print("Invalid memory settings.")
      memErr = 1
      break
    
    psi4_mol = psi4.core.Molecule.create_molecule_from_string(molString) 

    psi4_mol.update_geometry()
   
    # Given the test input, we are calculating the same calculations 10 times.
    # There should be a way to only calculate this once.
    #storing the single-point electronic energy in a variable
    try:
      ref_energy = psi4.energy("{}/{}".format(args.method, args.basis), 
          molecule=psi4_mol)
    except:
      print("Method does not exist.")
      methodErr = 1
      break

    # Store the one-body energy into molecule
    mol.setEnergy(ref_energy)
    print(mol.getEnergy())

    #Writing the one-body training set without 
    #parsing every output file in the end
    training_set_file.write(str(len(mol.atoms))+
        '\n'+str("%.8f" % ref_energy)+'\n'+molString+'\n')

training_set_file.close()

if memErr or methodErr or basisErr or fileErr:
    print("One or more errors have occurred. Exiting without saving.")
    sys.exit(1)

#Compression of the calculations directory
shutil.make_archive('calculations', 'zip', './', 'calculations')
