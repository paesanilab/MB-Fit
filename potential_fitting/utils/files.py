# external package imports
import os

def init_directory(directory_path):
    """
    Creates the directory (and any superdirectories) if they do not exist.

    Args:
        directory_path      - Local path to the directory.

    Returns:
        The same directory path as was passed in.
    """

    if not os.path.isdir(directory_path):
        os.makedirs(directory_path)

    return directory_path

def init_file(file_path):
    """
    Creates any directories that are needed to house the given file if they do not already exist.
    
    Args:
        file_path           - Local path to the file.

    Returns:
        The same file path as was passed in.
    """

    return os.path.join(init_directory(os.path.dirname(file_path)), os.path.basename(file_path))

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

    file_path = os.path.join(log_path, molecule.get_name(), "{}.{}".format(molecule.get_SHA1(), suffix))

    # make sure the required directories exist
    return init_file(file_path)
    
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
    return get_model_log_path(os.path.join(log_path, "optimizations"), method, basis, suffix)

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
    
    return get_model_log_path(os.path.join(log_path, "energy", "cp" if cp else "nocp"), method, basis, suffix)

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

    return get_model_log_path(os.path.join(log_path, "frequencies"), method, basis, suffix)
