# external package imports
import os
from glob import glob

# absolute module imports
from potential_fitting.molecule import Molecule
from potential_fitting.utils import SettingsReader
from potential_fitting.exceptions import ConfigMissingPropertyError

# local module imports
from .database import Database


def read_all_jobs(database_config_path, job_dir):
    """
    Searches the given directory for completed job directories and enters
    the results into the database.

    Any directory that starts with job_ will be considered a completed job directory.
    After data is entered into the database, these directories will be renamed from
    job_* to job_*_done.

    Args:
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        job_dir             - Local path the the directory to search.

    Returns:
        None.
    """
    calculation_results = []
    for directory in glob(job_dir + "/job_*"):
        if directory.endswith("done"):
            continue
        if not os.path.isdir(directory):
            continue
        calculation_results.append(read_job(directory + "/output.ini", directory + "/output.log"))

        if len(calculation_results) > 1000:
            with Database(database_config_path) as db:
                db.set_properties(calculation_results)
            calculation_results = []

    with Database(database_config_path) as db:
        db.set_properties(calculation_results)

    for directory in glob(job_dir + "/job_*"):
        if directory.endswith("done"):
            continue
        if not os.path.isdir(directory):
            continue

        i = 1

        job_dir = job_dir + "/job_{}_done".format(i)

        while os.path.exists(job_dir):
            i += 1

            job_dir = job_dir + "/job_{}_done".format(i)

        os.rename(directory, job_dir)


def read_job(job_dat_path, job_log_path):
    """
    Reads a completed job from its output file and log_file and enters the result into a database.
    
    Args:
        job_dat_path        - Local path to the .ini output file from this job.
        job_log_path        - Local path to the log file from this job.

    Returns:
        None.
    """

    data = SettingsReader(job_dat_path)

    xyz = data.get("molecule", "xyz")
    atom_counts = data.getlist("molecule", "atom_counts", int)
    charges = data.getlist("molecule", "charges", int)
    spins = data.getlist("molecule", "spins", int)
    symmetries = data.getlist("molecule", "symmetries", str)
    SMILES = data.get("molecule", "SMILES", str).split(",")
    names = data.getlist("molecule", "names", str)
    method = data.get("molecule", "method")
    basis = data.get("molecule", "basis")
    cp = data.get("molecule", "cp")
    use_cp = data.get("molecule", "use_cp")
    frag_indices = data.getlist("molecule", "frag_indices", int)

    molecule = Molecule.read_xyz(xyz, atom_counts, names, charges, spins, symmetries, SMILES)

    try:
        energy = data.getfloat("molecule", "energy")
        success = True
    except ConfigMissingPropertyError:
        energy = 0
        success = False

    with open(job_log_path, "r") as log_file:
        log_text = "\n".join(log_file.readlines())

    return molecule, method, basis, cp, use_cp, frag_indices, success, energy, log_text



