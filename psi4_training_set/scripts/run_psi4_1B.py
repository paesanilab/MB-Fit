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
# comment(mandatory, but can by empty)
# n atoms, and their xyz coordinates

# Check for initial passed-in command-line args
parser = argparse.ArgumentParser('Calculates single-body energy of molecules.')
parser.add_argument('memory', metavar='Memory', type=str, 
    default='500 MB', nargs='?', 
    help='Memory needed to make the calculation. Default is 500 MB.')
parser.add_argument('method', metavar='Method', type=str, default='HF',
    nargs='?', 
    help='The method to use for the calculation. Default is HF.')
parser.add_argument('basis', metavar='BasisSet', type=str, default='sto-3g',
    nargs='?', 
    help='The basis set to use for the calculation. Default is sto-3g.')
args = parser.parse_args()
print(args)

'''
if len(sys.argv) == 1:
    if len(sys.argv) != 1 and len(sys.argv) != 4:
    print("Usage: python run_psi4_1B.py (Memory Method BasisSet)")
    sys.exit(1)
'''

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
molecules = readfile('input.xyz')

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
      sys.exit(1)
    
    psi4_mol = psi4.core.Molecule.create_molecule_from_string(molString) 

    psi4_mol.update_geometry()
   
    # Given the test input, we are calculating the same calculations 10 times.
    # There should be a way to only calculate this once.
    #storing the single-point electronic energy in a variable
    ref_energy = psi4.energy('scf/3-21g', molecule=psi4_mol)

    # Store the one-body energy into molecule
    mol.setEnergy(ref_energy)

    #Writing the one-body training set without 
    #parsing every output file in the end
    training_set_file.write(str(len(mol.atoms))+
        '\n'+str("%.8f" % ref_energy)+'\n'+molString+'\n')

training_set_file.close()

#Compression of the calculations directory
shutil.make_archive('calculations', 'zip', './', 'calculations')
