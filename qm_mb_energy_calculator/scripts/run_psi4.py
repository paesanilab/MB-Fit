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

# TODO:
# Thread the calculation processes with psi4.set_num_thread

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
parser.add_argument('threads', type=str)

args = parser.parse_args()

training_set_file = open('./training_set.xyz', 'w')

poly_count = 0

# Read from input file
try:
    polymers = readfile(args.inputfile)
except:
    print("Invalid file. Exiting.")
    training_set_file.close()
    sys.exit(1)

# Perform a calculation for each polymer
for polymer in polymers:

    poly_count += 1
    final_energy = 0
    energy_str = ""

    # Deals with directory naming
    zeroes = int(math.ceil(len(str(poly_count)) / zfill_const) * zfill_const)
    dirname = str(poly_count).zfill(zeroes)

    # If a calculation does not exist, make a new file for it
    if not os.path.exists("{}{}/{}".format(pwd, calc_dir, dirname)):
        os.makedirs("{}{}/{}".format(pwd, calc_dir, dirname))

    # Designate the output file
    psi4.core.set_output_file("{}{}/{}{}".format(pwd, calc_dir, dirname, 
          output), False)
    
    # Perform calculations for all possible combinations in polymer
    for mol_comb in comb.make_combs(polymer.size()):

        # Determine whether this energy adds or subtracts.
        # For a trimer, trimer energy is added, dimer energy is subtracted,
        # and monomer energy is added.
        alternate = (-1)**(polymer.size() - len(mol_comb))

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
        '\n'+"%.8f"%final_energy + " " + energy_str+'\n'+str(polymer)+'\n')

training_set_file.close()

# Error message?
if mem_err or method_err or basis_err:
    print("One or more errors have occurred. Exiting without saving.")
    sys.exit(1)

#Compression of the calculations directory
shutil.make_archive('calculations', 'zip', pwd, calc_dir)
