# external package imports
import os

# absolute module imports
from potential_fitting.utils import SettingsReader, files
from potential_fitting.exceptions import ConfigMissingSectionError, ConfigMissingPropertyError

# local module imports
from .database import Database

def make_all_jobs(settings_path, client_name, job_dir, num_jobs):
    """
    Makes a Job file for each energy that still needs to be calculated in this Database.

    Args:
        settings_path       - Local path to the ".ini" file with relevent settings.
        client_name         - Name of the client that will perform these jobs
        job_dir             - Local path to the directory to place the job files in.
        num_jobs            - The number of jobs to generate.

    Returns:
        None.
    """

    # open the database
    with Database() as database:
        for molecule, method, basis, cp, frag_indices in database.get_all_calculations(client_name, calculations_to_do=num_jobs):
            write_job(settings_path, molecule, method, basis, cp, frag_indices, job_dir)

def write_job(settings_path, molecule, method, basis, cp, frag_indices, job_dir):
    """
    Makes a Job file for a specific Calculation.

    Args:
        settings_path       - Local path to the ".ini" file with relevent settings
        job                 - The Job object (see database.py) with the information needed to make a job
        job_dir             - Local path to the directory to place the job file in.

    Returns:
        None.
    
    """

    # parse settings file
    settings = SettingsReader(settings_path)

    i = 1

    file_path = job_dir + "/job_{}.py".format(i)

    while os.path.exists(file_path):
        i += 1

        file_path = job_dir + "/job_{}.py".format(i)

    files.init_file(file_path)

    with open(file_path, "w") as job_file, open(os.path.dirname(os.path.abspath(__file__)) + "/job_template.py", "r") as job_template:
        job_string = "".join(job_template.readlines())

        job_file.write(job_string.format(**{
            # TODO
            "whole_molecule": molecule.to_xyz().replace("\n", "\\n"),
            "charges": [frag.get_charge() for frag in molecule.get_fragments()],
            "spins": [frag.get_spin_multiplicity() for frag in molecule.get_fragments()],
            "symmetries": [frag.get_symmetry() for frag in molecule.get_fragments()],
            "atom_counts": [frag.get_num_atoms() for frag in molecule.get_fragments()],
            "names": [frag.get_name() for frag in molecule.get_fragments()],
            "total_atoms": molecule.get_num_atoms(),
            "molecule":     molecule.to_xyz(frag_indices, cp).replace("\n", "\\n"),
            "frag_indices": frag_indices,
            "method":       method,
            "basis":        basis,
            "cp":           cp,
            "num_threads":  settings.get("psi4", "num_threads"),
            "memory":       settings.get("psi4", "memory"),
            "format":       "{}"
        }))
