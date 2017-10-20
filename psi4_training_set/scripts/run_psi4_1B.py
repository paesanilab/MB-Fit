#!/usr/bin/env python3

import os, sys
import psi4
import pybel as pb
import shutil

from reader import readfile

#TODO
#Things to add later
#descriptive names for variables etc.
#reader for job control settings for psi4 like charge and multiplicity
#output storage and compression
#training set writing during/after calculation

# TODO
# Find a way to parse the test input file
# Cannot use BOTH pybel and psi4
# The format for xyz input file:
# amount of atoms in a molecule
# comment(mandatory, but can by empty)
# n atoms, and their xyz coordinates

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
    psi4.set_memory('500 MB')
    
    psi4_mol = psi4.core.Molecule.create_molecule_from_string(molString)  

    psi4_mol.update_geometry()
   
    #storing the single-point electronic energy in a variable
    ref_energy = psi4.energy('scf/3-21g', molecule=psi4_mol)

    #Writing the one-body training set without 
    #parsing every output file in the end
    training_set_file.write(str(len(mol.atoms))+
        '\n'+str("%.8f" % ref_energy)+'\n'+molString+'\n')

training_set_file.close()

#Compression of the calculations directory
shutil.make_archive('calculations', 'zip', './', 'calculations')
