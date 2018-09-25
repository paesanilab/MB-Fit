# absolute module imports
from potential_fitting.utils import SettingsReader
from potential_fitting.exceptions import ConfigMissingSectionError, ConfigMissingPropertyError

# local module imports
from .database import Database

def make_all_jobs(settings_path, database_path, job_dir):
    """
    Makes a Job file for each energy that still needs to be calculated in this Database.

    Args:
        settings_path       - Local path to the ".ini" file with relevent settings.
        database_path       - Local path to the database file. ".db" will be appended if it does not already end in
                ".db".
        job_dir             - Local path to the directory to place the job files in.

    Returns:
        None.
    """

    # open the database
    with Database(database_path) as database:

        for calculation in database.missing_energies():
            write_job(settings_path, calculation, job_dir)

def make_job(settings_path, database_path, job_dir):   
    """
    Makes a single Job file for an energy that still needs to be calculated in this Database.

    Args:
        settings_path       - Local path to the ".ini" file with relevent settings
        database_path       - Local path to the database file. ".db" will be appended if it does not already end in
                ".db".
        job_dir             - Local path to the directory to place the job file in.

    Returns:
        None.
    """

    # open the database
    with Database(database_path) as database:

        write_job(settings_path, database.get_missing_energy(), job_dir)

def write_job(settings_path, job, job_dir):
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

    with (open(job_dir + "/job_{}.py".format(job.job_id), "w") as job_file, open("job_template.py", "r") as
            job_template):
        job_string = "".join(job_template.readlines())

        job_file.write(job_string.format(**{
            "job_id":       job.job_id,
            "molecule":     job.molecule.to_xyz(job.fragments, job.cp).replace("\n", "\\n"),
            "method":       job.method,
            "basis":        job.basis,
            "num_threads":  settings.get("psi4", "num_threads"),
            "memory":       settings.get("psi4", "memory"),
            "format":       "{}"
        }))
