"""
A driver for the modules
"""

try:
    import tensorflow
except ImportError:
    pass

"""
LIBRARY IMPORTS
"""

#TODO: I don't really know what these imports do?
import os, sys

# Imports the configparser
import configparser
# Imports the json library
import json

"""
LOCAL MODULE IMPORTS
"""

# Calculates multi-body decomposition energies
import mbdecomp

# Holds information used in a molecule
import molecule

# Parses an xyz file into a molecule object
from molecule_parser import xyz_to_molecules

"""
CONFIG FILE INITIALIZATION
"""

# Config file, set by user, contains instructions for operation
config = configparser.ConfigParser(allow_no_value=False)

# read the conflig file into the configuration
config.read("settings.ini")


"""
FILEPATH INITIALIZATION
"""

# Filepath to the training set file
tr_set = './training_set.xyz'

# Filepath to the log file
log_name = './mbdecomp.log'

# Filepath to the input xyz file
input_path = config["driver"]["input"]

"""
OPEN FILES TO READ AND WRITE
"""

# open the input file
f = open(input_path, 'r')

# open the training set file
training_set_file = open(tr_set, 'w')

"""
INITIALIZE LOG FILE AND CHECK WHAT TYPE OF LOGGING IS TO BE DONE
"""

# True if a log file is to be written
write_log = False

# Open the log file if cluster_anaylsis
if config["MBdecomp"].getboolean("cluster_analysis"):
            
    # Check if the append flag is set to true
    if config["MBdecomp"].getboolean("append"):
        # If both append and overwrite flags are True, print a warning message
        if config["MBdecomp"].getboolean("overwrite"):
            print("Warning: Both flags selected. Will append " + 
                "to file if it exists.")
        # Open the log file in append mode
        log = open(log_name, 'a')
        write_log = True
    # Check if overwrite flag is set to true
    elif config["MBdecomp"].getboolean("overwrite"):
        # Open the log file in write mode
        log = open(log_name, 'w')
        write_log = True

    # If file exists and both flags are set to false, no not log
    elif os.path.isfile(log_name):
        print("Error: File exists and don't know how to deal with it."+
            " No log file will be produced.")
        
'''
Notes for output file:
Let user decide what to do with log file.
By default, if a log file already exists, terminate.
We can also let the user append or overwrite the existing file.
'''

"""
RUN CALCULATIONS ON THE USER'S INPUT
"""

# parse array of molecules from input xyz file
molecules = xyz_to_molecules(f)

# for each molecule from the input xyz file...
for molecule in molecules:
    # calculate energy
    energy = mbdecomp.get_nmer_energies(molecule, config)
    # calculate mb_energies
    molecule.mb_energies = mbdecomp.mbdecomp(molecule.nmer_energies)
    # write output to training set file
    training_set_file.write(str(molecule.get_num_atoms()) + '\n' + str(energy) + '\n' + molecule.to_xyz() + '\n')

    # if log is enabled, write the log
    if write_log:
        json_output = {}
        
        # write molecule's xyz format to log file
        log.write(molecule.to_xyz() + '\n')
        # write molecule's fragment energies to log file
        log.write(molecule.log_frag_energy())
        log.write("N-body energies:\n")
        # write molecule's multibody energy to log file
        log.write(molecule.log_mb_energy(config["MBdecomp"].getint("max_nbody_energy")))
        
        if config["MBdecomp"].getboolean("kbody_energy"):
            log.write("K-body energies:\n")
            # get kbody energies
            k_output = mbdecomp.get_kbody_energies(molecule)
            k_dict = k_output[0]
            k_diff = k_output[1]
            for index in k_dict.keys():
                log.write("{}: {}\n".format(index, "%.8f"%k_dict[index]))
            log.write("K-body differences:\n")
            for bodies in range(len(k_diff)):
                log.write("V_{}B - K_{}B: {}\n".format(bodies+1, bodies+1, "%.8f"%k_diff[bodies]))
        log.write("--------------\n")
        
        # Build json dictionary
        
        # put molecule's xyz format into json
        json_output["molecule"] = molecule.to_xyz()
        
        # put fragment energies into json
        json_output["frag_energies"] = {"E{}".format(key): molecule.energies[key] for key in molecule.energies.keys()}
        
        # put nbody energies into json
        json_output["n-body_energies"] = {"V_{}B".format(molecule.mb_energies.index(energy) + 1): energy for energy in molecule.mb_energies}
        
        # put kbody energies into json
        json_output["k-body_energies"] = k_dict
        with open("json_output.json", 'w') as json_file:
            json.dump(json_output, json_file, indent=4)


            
        

"""
OLD CODE, COMMENTED
while while_read:
    if while_read:
        count += 1
        #print(mol_from_xyz)
        #print(count)
        energy = mbdecomp.get_nmer_energies(mol_from_xyz, config)
        mol_from_xyz.mb_energies = mbdecomp.mbdecomp(
            mol_from_xyz.nmer_energies[::-1])
        training_set_file.write(str(mol_from_xyz.natoms)+
            '\n'+str(energy)+'\n'+str(mol_from_xyz)+'\n')
        
        if write_log:
            # Create a dictionary (for JSON) that collects all output info
            # How to check the uniqueness of a molecule?
            # One proposal: generate and check MD5 checksum
            json_output = {}
            
            # Continue to write this in a log file
            log.write(str(mol_from_xyz)+'\n')
            log.write(mol_from_xyz.write_frag_energy())
            log.write("N-body energies:\n")
            log.write(mol_from_xyz.write_mb_energy(
                  config["MBdecomp"].getint("max_nbody_energy")))
            if config["MBdecomp"].getboolean("kbody_energy"):
                  log.write("K-body energies:\n")
                  k_output = mbdecomp.get_kbody_energies(mol_from_xyz)
                  k_dict = k_output[0]
                  k_diff = k_output[1]
                  for index in k_dict.keys():
                      log.write("{}: {}\n".format(index, "%.8f"%k_dict[index]))
                  log.write("K-body differences:\n")
                  for bodies in range(len(k_diff)):
                      log.write("V_{}B - K_{}B: {}\n".format(bodies+1,
                          bodies+1, "%.8f"%k_diff[bodies]))
            log.write("--------------\n")

            # Place some entries of output into dictionary
            json_output["molecule"] = str(mol_from_xyz)
            # nmer_energies
            json_output["frag_energies"] = {"E{}".format(key): 
                mol_from_xyz.energies[key] for key in 
                mol_from_xyz.energies.keys()}        
            json_output["n-body_energies"] = {
                "V_{}B".format(mol_from_xyz.mb_energies.index(energy)+1):
                energy for energy in mol_from_xyz.mb_energies}
            json_output["k-body_energies"] = k_dict
            with open('json_output.json', 'w') as fp:
                json.dump(json_output, fp, indent=4)
"""

# We make an exception here in case there is no log file
try:
    log.close()
except:
    pass

f.close()

training_set_file.close()
