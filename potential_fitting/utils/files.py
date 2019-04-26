# external package imports
import os
import datetime
from enum import Enum

# absolute package imports
#from potential_fitting.exceptions import FileExistsError, InvalidValueError
from potential_fitting.exceptions import (LibraryNotAvailableError, LibraryCallError, NoSuchLibraryError,
        ConfigMissingSectionError, ConfigMissingPropertyError, CommandExecutionError)

def init_directory(directory_path):
    """
    Creates the directory (and any superdirectories) if they do not exist.

    Args:
        directory_path      - Local path to the directory.

    Returns:
        The same directory path as was passed in.
    """
    if directory_path == "":
        directory_path = "."
    if not os.path.isdir(directory_path):
        os.makedirs(directory_path)

    return directory_path

class OverwriteMethod(Enum):
    OVERWRITE = 1
    BACKUP = 2
    CRASH = 3
    NONE = 4 # special overwrite method for when desired behavior is modification of an existing file or silent
             # silent overwritting of the file.

    def get_from_settings(settings):
        overwrite_method = settings.get("files", "overwrite_method", "backup")
        if overwrite_method.lower() == "overwrite":
            return OverwriteMethod.OVERWRITE
        if overwrite_method.lower() == "backup":
            return OverwriteMethod.BACKUP
        if overwrite_method.lower() == "crash":
            return OverwriteMethod.CRASH
        if overwrite_method.lower() == "none":
            raise InvalidValueError("overwrite_method", overwrite_method, "none is a special overwrite method and "
                    + "should not be used by the user.")

        raise InvalidValueError("overwrite_method", overwrite_method, "overwrite_method must be a valid "
                + "OverwriteMethod. See files.OverwriteMethod")
        


def init_file(file_path, overwrite_method = OverwriteMethod.BACKUP):
    """
    Creates any directories that are needed to house the given file if they do not already exist.
    
    Args:
        file_path           - Local path to the file.

    Returns:
        The same file path as was passed in.
    """

    file_path = os.path.join(init_directory(os.path.dirname(file_path)), os.path.basename(file_path))

    if os.path.isfile(file_path):
        if overwrite_method == OverwriteMethod.OVERWRITE:
            print("File {} already exists, but overwrite_method is set to OVERWRITE, so it is being overwritten.".format(file_path))

        elif overwrite_method == OverwriteMethod.BACKUP:
            i = 1

            backup_path = file_path + ".backup-{}".format(i)
            while os.path.isfile(backup_path):
                i += 1
                backup_path = file_path + ".backup-{}".format(i)

            print("File {0} already exists, moving existing {0} to {1} to make way for new file.".format(file_path,
                    backup_path))

            backup_path = init_file(backup_path)

            os.rename(file_path, backup_path)

        elif overwrite_method == OverwriteMethod.CRASH:
            raise FileExistsError(file_path)

        elif overwrite_method == OverwriteMethod.NONE:
            pass

        else:
            raise InvalidValueError("overwrite_method", overwrite_method, "must be a valid OverwriteMethod. "
                    + "See files.OverwriteMethod")
    return file_path

def get_molecule_log_path(log_path, molecule, suffix):
    """
    Returns the path to a log file for this molecule and creates any directories needed to house the log file.

    Args:
        log_path            - Local path to the logging directory to make this log file in.
        molecule            - The molecule to get the log file of.
        suffix              - The type of log file, should be one if "in", "log", or "out".

    Returns:
        The log file for the given molecule with the given suffix.
    """

    if not "." in suffix:
        file_path = os.path.join(log_path, molecule.get_name(), "{}.{}".format(molecule.get_SHA1()[-8:], suffix))
    else:
        file_path = os.path.join(log_path, molecule.get_name(), "{}{}".format(molecule.get_SHA1()[-8:], suffix))

    # make sure the required directories exist
    return init_file(file_path, overwrite_method = OverwriteMethod.NONE)
    
def get_model_log_path(log_path, molecule, method, basis, suffix):
    """
    Returns the path to a log file for this molecule, method, and basis.

    Args:
        log_path            - Local path to the logging directory to make this log file in.
        molecule            - The molecule to get the log file of.
        method              - The method of this log path.
        basis               - The basis of this log path.
        suffix              - The type of log file, should be one if "in", "log", or "out".

    Returns:
        The log file for the given molecule, method, and basis with the given suffix
    """
    
    return get_molecule_log_path(os.path.join(log_path, method, basis), molecule, suffix)

def get_optimization_log_path(log_path, molecule, method, basis, suffix):
    """
    Returns the path to an optimization log file for this molecule.

    Args:
        log_path            - Local path to the logging directory to make this log file in.
        molecule            - The molecule to get the log file of.
        method              - The method of this log path.
        basis               - The basis of this log path.
        suffix              - The type of log file, should be one if "in", "log", or "out".

    Returns:
        The log file for the optimization of this molecule using the given method/basis with the given suffix.
    """
    return get_model_log_path(os.path.join(log_path, "optimizations"), molecule, method, basis, suffix)

def get_energy_log_path(log_path, molecule, method, basis, cp, suffix):
    """
    Returns the path to an energy log file for this molecule.

    Args:
        log_path            - Local path to the logging directory to make this log file in.
        molecule            - The molecule to get the log file of.
        method              - The method of this log path.
        basis               - The basis of this log path.
        cp                  - The counterpoise correction of this log path.
        suffix              - The type of log file, should be one if "in", "log", or "out".

    Returns:
        The log file for the energy calculation of this molecule using the given method/basis/cp with the given suffix.
    """
    
    return get_model_log_path(os.path.join(log_path, "energy", "cp" if cp else "nocp"), molecule, method, basis,
            suffix)

def get_frequencies_log_path(log_path, molecule, method, basis, suffix):
    """
    Returns the path to a frequencies log file for this molecule.

    Args:
        log_path            - Local path to the logging directory to make this log file in.
        molecule            - The molecule to get the log file of.
        method              - The method of this log path.
        basis               - The basis of this log path.
        suffix              - The type of log file, should be one if "in", "log", or "out".

    Returns:
        The log file for the frequencies calculation of this molecule using the given method/basis with the given
        suffix.
    """

    return get_model_log_path(os.path.join(log_path, "frequencies"), molecule, method, basis, suffix)

def get_qchem_input_string(molecule, fragment_indicies, model, cp, settings):
    method, basis = model.split("/")
    # initialize qchem input string
    qchem_input = "";

    # molecule format
    qchem_input += "$molecule\n"

    # charge and spin multiplicity
    qchem_input += "{} {}\n".format(molecule.get_charge(fragment_indicies),
            molecule.get_spin_multiplicity(fragment_indicies))

    # atoms in the molecule
    # might need to add whitespace before each line?
    qchem_input += molecule.to_xyz(fragment_indicies, cp) + "\n"

    qchem_input += "$end\n"

    # Q-chem settings
    qchem_input += "$rem\n"
    qchem_input += "jobtype " + "sp" + "\n"
    qchem_input += "method " + model.split('/')[0] + "\n"
    qchem_input += "basis " + model.split('/')[1] + "\n"

    try:
        qchem_input += "ecp " + settings.get("qchem", "ecp") + "\n"
    except (ConfigMissingSectionError, ConfigMissingPropertyError):
        pass

    qchem_input += "$end"

    return qchem_input

def write_qchem_input(settings, qchem_in_path, jobtype, template_input = ""):
    with open(qchem_in_path, "w") as qchem_in_file:

        # indicate to qchem that the molecule begins here
        qchem_in_file.write("$molecule\n")

        # tells qchem the charge and spin multiplicity of this molecule
        qchem_in_file.write("{} {}\n".format(molecule.get_charge(), molecule.get_spin_multiplicity()))
        
        
        # add molecule's xyz to qchem input
        qchem_in_file.write(molecule.to_xyz(fragment_indicies, cp) + "\n")

        # indicate to qchem that the molecule ends here
        qchem_in_file.write("$end\n")

        # indicate to qchem that settings start here
        qchem_in_file.write("$rem\n")
        

        # indicate what type of calculation to run
        qchem_in_file.write("jobtype {}\n".format(jobtype))

        # indicate the method and basis to use in the calculation
        qchem_in_file.write("method {}\n".format(model.split('/')[0]))
        qchem_in_file.write("basis {}\n".format(model.split('/')[1]))

        try:
            # if the user whishes to specify ecp, then specify an ecp
            qchem_in_file.write("ecp {}\n".format(settings.get("qchem", "ecp")))

        except (ConfigMissingSectionError, ConfigMissingPropertyError):
            # if the user did not specify ecp, then ommit it
            pass

        # if the user included a special input template with additional arguments, include it.
        qchem_in_file.write(template_input)

        # indicate to qchem that settings end here
        qchem_in_file.write("$end")
