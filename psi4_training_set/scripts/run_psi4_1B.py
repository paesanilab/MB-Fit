#!/usr/bin/env python3

import os, sys
import psi4
import pybel as pb
import shutil
import argparse
import math

from reader import readfile
import calculate

zfill_const = 4
pwd = './'
calcDir = 'calculations'
output = '/output.dat'

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
memErr = methodErr = basisErr = 0

'''
Subclass and function to override argparse
'''
class optionParse(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_lines):
        return arg_lines.split()[1:]

# Check for initial passed-in command-line args
parser = optionParse(fromfile_prefix_chars='@')

# Add arguments and parse
parser.add_argument('memory', type=str)
parser.add_argument('method', type=str) 
parser.add_argument('basis', type=str)
parser.add_argument('charge', type=int)
parser.add_argument('spin', type=int)
parser.add_argument('convergence', type=int)
parser.add_argument('cycles', type=int)
parser.add_argument('threshold', type=int)
parser.add_argument('inputfile', type=str)
args = parser.parse_args()

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
    # Zeroes should no longer be hardcoded - Derek

# Temporary solution to pybel conflict. See reader.py
try:
    molecules = readfile(args.inputfile)
except:
    print("Invalid file. Exiting.")
    training_set_file.close()
    sys.exit(1)

# Perform a calculation for each molecule
for mol in molecules:
    molCount += 1
    molString = mol.toString()

    zeroes = int(math.ceil(len(str(molCount)) / zfill_const) * zfill_const)
    dirname = str(molCount).zfill(zeroes)

    # If a calculation does not exist, make a new file for it
    if not os.path.exists("{}{}{}".format(pwd, calcDir, dirname)):
        os.makedirs("{}{}{}".format(pwd, calcDir, dirname))

    # Designate the output file
    psi4.core.set_output_file("{}{}{}/{}".format(pwd, calcDir, dirname, 
          output), False)
    
    # Redirect calculation
    ref_energy = calculate.calculate(mol, args)

    #Writing the one-body training set without 
    #parsing every output file in the end
    training_set_file.write(str(len(mol.atoms))+
        '\n'+str("%.8f" % ref_energy)+'\n'+molString+'\n')

training_set_file.close()

if memErr or methodErr or basisErr:
    print("One or more errors have occurred. Exiting without saving.")
    sys.exit(1)

#Compression of the calculations directory
shutil.make_archive('calculations', 'zip', './', 'calculations')
