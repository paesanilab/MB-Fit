from subprocess import call
import json
from mbfit import molecule
from mbfit.calculator import mbdecomp
from mbfit.molecule import xyz_to_molecules

from mbfit.utils import SettingsReader
from mbfit.exceptions import InvalidValueError, ParsingError, XYZFormatError, InconsistentValueError

# Config file, set by user, contains instructions for operation
settings = SettingsReader("driver_settings.ini")

# Filepath to the training set file
tr_set = './training_set.xyz'

# Filepath to the log file
log_name = './mbdecomp.log'

# Filepath to the input xyz file
input_path = settings.get("driver", "input")

"""
INITIALIZE LOG FILE AND CHECK WHAT TYPE OF LOGGING IS TO BE DONE
"""

# True if a log file is to be written
write_log = False

# Open the log file if cluster_anaylsis
if settings.getboolean("MBdecomp", "cluster_analysis"):
    
    # Check if the log type is overwrite
    if settings.get("driver", "log_type") == "overwrite":
        write_log = True
        # clear the log file if set to write
        with open(log_name, "w") as log_file:
            pass
    # Check if the log type is append
    elif settings.get("driver", "log_type") == "append":
        write_log = True
    # Check if log type is invalid
    elif settings.get("driver", "log_type") != "none":
        raise InvalidValueError("[driver][log_type]", settings.get("driver", "log_type"), "one of 'overwrite', 'append' or 'none'")

"""
RUN CALCULATIONS ON THE USER'S INPUT
"""

# parse array of molecules from input xyz file
try:
    molecules = xyz_to_molecules(input_path, settings)
except (XYZFormatError, InconsistentValueError) as e:
    raise ParsingError(input_path, str(e)) from None

with open(tr_set, "w") as training_set_file:
    # for each molecule from the input xyz file...
    for molecule in molecules:

        mol_nfrags = molecule.get_num_fragments()

        # calculate energy
        energy = mbdecomp.get_nmer_energies("driver_settings.ini", molecule)
        
        # calculate mb_energies
        molecule.mb_energies = mbdecomp.mbdecomp(molecule.nmer_energies)

        # write output to training set file
        training_set_file.write(str(molecule.get_num_atoms()) + "\n" +
            str(energy) + "\n" + molecule.to_xyz() + "\n")

        # if log is enabled, write the log
        if write_log:
            with open(log_name, "a") as log:
                json_output = {}
                
                # write molecule's xyz format to log file
                log.write("Molecule: \n" + molecule.to_xyz() + "\n")
                # write molecule's fragment energies to log file
                log.write("Fragment Energies: \n" + molecule.log_frag_energy())
                log.write("N-body energies:\n")
                # write molecule's multibody energy to log file
                log.write(molecule.log_mb_energy(settings.getint("MBdecomp", "max_nbody_energy")))
                
                if settings.getboolean("MBdecomp", "kbody_energy"):
                    log.write("Individual N-body energies:\n")
                    # get kbody energies
                    k_output = mbdecomp.get_kbody_energies(molecule)

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
                json_output["model"] = settings.get("model", "method") + "/" + settings.get("model", "basis")
                
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
