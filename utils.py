import os

import subprocess

from exceptions import CommandNotFoundError, CommandExecutionError

def init_directory(directory_path):
    """
    Creates the directory (and any superdirectories) if they do not exist

    Args:
        directory_path - the path to the directory

    Returns:
        the same directory path as was passed in
    """

    if not os.path.isdir(directory_path):
        os.makedirs(directory_path)

    return directory_path

def init_file(file_path):
    """
    Creates any directories that are needed to house the given file if they do not already exist
    
    Args:
        file_path - the path to the file

    Returns:
        the same file path as was passed in
    """

    return os.path.join(init_directory(os.path.dirname(file_path)), os.path.basename(file_path))

def get_molecule_log_path(log_path, molecule, suffix):
    """
    Returns the path to a log file for this molecule

    Args:
        log_path    - the path to the logging directory to make this log file in
        molecule    - the molecule to get the log file of
        suffix      - the type of log file, should be one if "in", "log", or "out"

    Returns:
        the log file for the given molecule
    """

    file_path = os.path.join(log_path, molecule.get_name(), "{}.{}".format(molecule.get_SHA1(), suffix))

    # make sure the required directories exist
    init_file(path)
    
    return file_path

def get_model_log_path(log_path, molecule, method, basis, cp, suffix):
    """
    Returns the path to a log file for this molecule, method, basis, and cp

    Args:
        log_path    - the path to the logging directory to make this log file in
        molecule    - the molecule to get the log file of
        method      - the method of this log path
        basis       - the basis of this log path
        cp          - the counterpoise correction of this log path
        suffix      - the type of log file, should be one if "in", "log", or "out"

    Returns:
        the log file for the given molecule, method, basis, cp suffix
    """
    
    return get_molecule_log_path(os.path.join(log_path, method, basis, cp), molecule, suffix)

def get_optimization_log_path(log_path, molecule, method, basis, cp, suffix):
    """
    Returns the path to an optimization log file for this molecule

    Args:
        log_path    - the path to the logging directory to make this log file in
        molecule    - the molecule to get the log file of
        method      - the method of this log path
        basis       - the basis of this log path
        cp          - the counterpoise correction of this log path
        suffix      - the type of log file, should be one if "in", "log", or "out"

    Returns:
        the log file for the optimization of this molecule
    """
    return get_model_log_path(os.path.join(log_path, "optimizations"), method, basis, cp, suffix)

def get_energy_log_path(log_path, molecule, method, basis, cp, suffix):
    """
    Returns the path to an energy log file for this molecule

    Args:
        log_path    - the path to the logging directory to make this log file in
        molecule    - the molecule to get the log file of
        method      - the method of this log path
        basis       - the basis of this log path
        cp          - the counterpoise correction of this log path
        suffix      - the type of log file, should be one if "in", "log", or "out"

    Returns:
        the log file for the energy calculation of this molecule
    """
    
    return get_model_log_path(os.path.join(log_path, "energy"), method, basis, cp, suffix)

def get_frequencies_log_path(log_path, molecule, method, basis, cp, suffix):
    """
    Returns the path to a frequencies log file for this molecule

    Args:
        log_path    - the path to the logging directory to make this log file in
        molecule    - the molecule to get the log file of
        method      - the method of this log path
        basis       - the basis of this log path
        cp          - the counterpoise correction of this log path
        suffix      - the type of log file, should be one if "in", "log", or "out"

    Returns:
        the log file for the frequencies calculation of this molecule
    """

    return get_model_log_path(os.path.join(log_path, "frequencies"), method, basis, cp, suffix)

def get_1b_config_log_path(log_path, molecule, suffix):
    """
    Returns the path to the 1b config log file for this molecule

    Args:
        log_path    - the path to the logging directory to make this log file in
        molecule    - the molecule to get the log file of
        method      - the method of this log path
        basis       - the basis of this log path
        cp          - the counterpoise correction of this log path
        suffix      - the type of log file, should be one if "in", "log", or "out"

    Returns:
        the log file for the 1b config log file of this molecule
    """

    return get_molecule_log_path(os.path.join(log_path, "configuration_generator"), molecule, suffix)

def sys_call(command, *args):
    """
    Performs a system call with the given command and arguments

    A CommandNotFoundError will be raised if the command is not found (exit code 127)

    A CommandExecutionError will be raised if the exit code is not 0

    Args:
        command - the command to be executed
        args    - any arguments to follow the command

    Returns:
        None
    """

    try:
        subprocess.run([command, *args], stderr = subprocess.PIPE, check = True)
    except subprocess.CalledProcessError as e:
        if e.returncode == 127:
            raise CommandNotFoundError(command)
        raise CommandExecutionError(command, e.cmd, e.returncode, e.stderr)
    except FileNotFoundError:
        raise CommandNotFoundError(command)
