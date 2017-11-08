#!/usr/bin/env python3

import os, sys
import psi4
import pybel as pb
import shutil
import argparse
import math

from reader import readfile
import calculate
import comb

zfill_const = 4
pwd = './'
calc_dir = 'calculations'
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
mem_err = method_err = basis_err = 0

class OptionParse(argparse.ArgumentParser):
    """ Subclass and function to override default argparse parser """
    def convert_arg_line_to_args(self, arg_lines):
        return arg_lines.split()[1:]

# Check for initial passed-in command-line args
parser = OptionParse(fromfile_prefix_chars='@')

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


poly_count = 0

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
'''
try:
    polymers = readfile(args.inputfile)
except:
    print("Invalid file. Exiting.")
    training_set_file.close()
    sys.exit(1)
'''
polymers = readfile(args.inputfile)
# print(polymers)

# Perform a calculation for each molecule
for polymer in polymers:

    poly_count += 1
    final_energy = 0
    energy_str = ""

    zeroes = int(math.ceil(len(str(poly_count)) / zfill_const) * zfill_const)
    dirname = str(poly_count).zfill(zeroes)

    # If a calculation does not exist, make a new file for it
    if not os.path.exists("{}{}/{}".format(pwd, calc_dir, dirname)):
        os.makedirs("{}{}/{}".format(pwd, calc_dir, dirname))

    # Designate the output file
    psi4.core.set_output_file("{}{}/{}{}".format(pwd, calc_dir, dirname, 
          output), False)
    
    for mol_comb in comb.make_combs(polymer.size()):

        # Determine whether this energy adds or subtracts.
        # For a trimer, trimer energy is added, dimer energy is subtracted,
        # and monomer energy is added.
        #print(polymer.size())
        #print(len(mol_comb))
        alternate = (-1)**(polymer.size() - len(mol_comb))
        #print(alternate)

        # Redirect calculation
        ref_energy = calculate.calculate(polymer, mol_comb, args)
        
        # Add to final energy
        ref_energy *= alternate
        print("Energy for this calculation: {}".format(ref_energy))
        print("Number of molecules: {}".format(len(mol_comb)))
        energy_str += "%.8f"%ref_energy + " "
        final_energy += ref_energy

    # Print out the final energy calculation based on several calculations
    print("Final energy: {}".format(final_energy))

    #Writing the one-body training set without 
    #parsing every output file in the end
    training_set_file.write(str(polymer.total_atoms())+
        '\n'+energy_str+'\n'+str(polymer)+'\n')

training_set_file.close()

if mem_err or method_err or basis_err:
    print("One or more errors have occurred. Exiting without saving.")
    sys.exit(1)

#Compression of the calculations directory
shutil.make_archive('calculations', 'zip', pwd, calc_dir)
