# external package imports
import os
from glob import glob

# absolute module imports
from potential_fitting.molecule import Molecule
from potential_fitting.utils import SettingsReader
from potential_fitting.exceptions import ConfigMissingPropertyError

# local module imports
from .database import Database

def read_all_jobs(job_dir):
    calculation_results = []
    for directory in glob(job_dir + "/job_*"):
        print(directory)
        if directory.endswith("done"):
            continue
        calculation_results.append(read_job(directory + "/output.ini", directory + "/output.log"))

        i = 1

        job_dir = "job_{}_done".format(i)

        while os.path.exists(job_dir):
            i += 1

            job_dir = "job_{}_done".format(i)

        os.rename(directory, job_dir)

        if len(calculation_results) > 1000:
            with Database() as db:
                db.set_properties(calculation_results)

    with Database() as db:
        db.set_properties(calculation_results)


def read_job(job_dat_path, job_log_path):
    """
    Reads a completed job from its output file and enters the result into a database.
    
    Args:
        database_path       - Local path to the file where the database is stored. ".db" will be appended if it does
                not already end in "db".
        job_path            - Local path to the job_<id>.out output file to enter into the datbase.
        job_log_path        - Local path to the log file from this job.

    Returns:
        None
    """

    data = SettingsReader(job_dat_path)

    xyz = data.get("molecule", "xyz")
    atom_counts = data.getlist("molecule", "atom_counts", int)
    charges = data.getlist("molecule", "charges", int)
    spins = data.getlist("molecule", "spins", int)
    symmetries = data.getlist("molecule", "symmetries", str)
    names = data.getlist("molecule", "names", str)
    method = data.get("molecule", "method")
    basis = data.get("molecule", "basis")
    cp = data.get("molecule", "cp")
    frag_indices = data.getlist("molecule", "frag_indices", int)

    print(symmetries)

    molecule = Molecule().read_xyz(xyz, atom_counts, names, charges, spins, symmetries)

    try:
        energy = data.getfloat("molecule", "energy")
        success = True
    except (ConfigMissingPropertyError):
        energy = 0
        success = False


    log_text = ""
    with open(job_log_path, "r") as log_file:
        log_text = "\n".join(log_file.readlines())


    return molecule, method, basis, cp, frag_indices, success, energy, log_text



