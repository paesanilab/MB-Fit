"""
A driver for the modules
"""

"""
LIBRARY IMPORTS
"""

# System-related imports for all your file needs
import os, sys
from subprocess import call

# Data compression modules
import pickle

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
config.read("driver_settings.ini")


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
    
    # Check if the log type is overwrite
    if config["driver"]["log_type"] == "overwrite":
        log = open(log_name, 'w')
        write_log = True
    # Check if the log type is append
    elif config["driver"]["log_type"] == "append":
        # Open the log file in write mode
        log = open(log_name, 'w')
        write_log = True
    # Check if log type is invalid
    elif config["driver"]["log_type"] != "none":
        raise ValueError("Ivalid log_type {}".format(config["driver"]["log_type"]))

"""
RUN CALCULATIONS ON THE USER'S INPUT
"""

# parse array of molecules from input xyz file
molecules = xyz_to_molecules(f, config)

# for each molecule from the input xyz file...
for molecule in molecules:

    mol_nfrags = molecule.get_num_fragments()

    # calculate energy
    energy = mbdecomp.get_nmer_energies(molecule, config)
    print(molecule.nmer_energies)
    
    # calculate mb_energies
    molecule.mb_energies = mbdecomp.mbdecomp(molecule.nmer_energies)

    # write output to training set file
    training_set_file.write(str(molecule.get_num_atoms()) + "\n" +
        str(energy) + "\n" + molecule.to_xyz() + "\n")

    # if log is enabled, write the log
    if write_log:
        json_output = {}
        
        # write molecule's xyz format to log file
        log.write("Molecule: \n" + molecule.to_xyz() + "\n")
        # write molecule's fragment energies to log file
        log.write("Fragment Energies: \n" + molecule.log_frag_energy())
        log.write("N-body energies:\n")
        # write molecule's multibody energy to log file
        log.write(molecule.log_mb_energy(config["MBdecomp"].getint("max_nbody_energy")))
        
        if config["MBdecomp"].getboolean("kbody_energy"):
            log.write("Individual N-body energies:\n")
            # get kbody energies
            k_output = mbdecomp.get_kbody_energies(molecule)

            #TODO: understand and comment this            

            k_dict = k_output[0]
            k_diff = k_output[1]
            for index in k_dict.keys():
                log.write("{}: {}\n".format(index, "%.8f"%k_dict[index]))
            log.write("N-body differences:\n")
            for bodies in range(len(k_diff)):
                log.write("V_{}B - K_{}B: {}\n".format(bodies+1, bodies+1, "%.8f"%k_diff[bodies]))
        log.write("--------------\n")
        
        # Build json dictionary

        json_output["ID"] = molecule.get_SHA1()

        # put molecule's xyz format into json
        json_output["molecule"] = molecule.to_xyz()

        # put the computational model into json
        json_output["model"] = config["model"]["method"]+"/"+config["model"]["basis"]
        
        # put fragment energies into json
        json_output["frag_energies"] = {"E{}".format(key): molecule.energies[key] for key in molecule.energies.keys()}
        
        # put nbody energies into json
        json_output["n-body_energies"] = {"V_{}B".format(molecule.mb_energies.index(energy) + 1): energy for energy in molecule.mb_energies}
        
        # put kbody energies into json
        json_output["individual_n-body_energies"] = k_dict
        with open("json_output.json", 'w') as json_file:
            json.dump(json_output, json_file, indent=4)
        
# We make an exception here in case there is no log file
if write_log:
    log.close()

f.close()

training_set_file.close()
