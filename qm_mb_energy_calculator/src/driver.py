"""
A driver for the modules
"""

#For some reason this import gives a warning.
try:
    import tensorflow
except ImportError:
    pass

import os, sys
import configparser
import json

import mbdecomp
import molecule

# 
tr_set = './training_set.xyz'
# 
log_name = './mbdecomp.log'
# 
overwrite = 'w'
# 
append = 'a'
# Whether or not to write a log file
write_log = False

# config file, set by user, contains instructions for operation
config = configparser.ConfigParser(allow_no_value=False)
# read the conflig file into the configuration
config.read("settings.ini")

# open the input file
f = open(config["driver"]["input"], 'r')
# open the training set file
training_set_file = open('./training_set.xyz', 'w')


# Open the log file here
if config["MBdecomp"].getboolean("cluster_analysis"):
            
    # Check if the file already exists, and do something about it
    if config["MBdecomp"].getboolean("append"):
        if config["MBdecomp"].getboolean("overwrite"):
            print("Warning: Both flags selected. Will append " + 
                "to file if it exists.")
        log = open(log_name, 'a')
        write_log = True
    elif config["MBdecomp"].getboolean("overwrite"):
        log = open(log_name, 'w')
        write_log = True

    # If file exists and both flags are set to false
    elif os.path.isfile(log_name):
        print("Error: File exists and don't know how to deal with it."+
            " No log file will be produced.")
        


while_read = 1
count = 0
'''
Notes for output file:
Let user decide what to do with log file.
By default, if a log file already exists, terminate.
We can also let the user append or overwrite the existing file.
'''
while while_read:
    mol_from_xyz = molecule.Molecule()
    while_read = mol_from_xyz.read_xyz(f)
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


# We make an exception here in case there is no log file
try:
    log.close()
except:
    pass

f.close()

training_set_file.close()
