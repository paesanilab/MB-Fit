from database import Database
from potential_fitting.utils import SettingsReader
from potential_fitting.exceptions import ConfigMissingSectionError, ConfigMissingPropertyError

def make_all_jobs(settings_file, database_name, job_dir):
    """
    Makes all the jobs that still need to be performed in this database

    Args:
        settings_file - .ini file with relevent settings
        database_name - file path to where the database is stored
        job_dir - directory to place the jobs in

    Returns:
        None
    """
    
    # open the database
    database = Database(database_name)
    with Database(database_name) as database:

        for calculation in database.missing_energies():
            write_job(settings_file, calculation, job_dir)

def make_job(settings_file, database_name, job_dir):   
    """
    Makes a single job that still needs to be performed in this database

    Args:
        settings_file - .ini file with relevent settings
        database_name - file path to where the database is stored
        job_dir - directory to place the job in

    Returns:
        None
    """
    
    # add .db to the database name if it doesn't already end in .db
    if database_name[-3:] != ".db":
        database_name += ".db"

    # open the database
    with Database(database_name) as database:

        write_job(settings_file, database.get_missing_energy(), job_dir)

def write_job(settings_file, calculation, job_dir):
    """
    Makes a job from a specific Calculation

    file will be located at job_dir/job_<id>.py
    
    Args:
        settings_file - .ini file with relevent settings
        calculation - the Calculation object (see database.py) with the information needed to make a job
        job_dir - directory to place the job in

    Returns:
        None
    
    """
    # parse settings file
    settings = SettingsReader(settings_file)

    with open(job_dir + "/job_{}.py".format(calculation.job_id), "w") as job_file, open("job_template.py", "r") as job_template:
        job_string = "".join(job_template.readlines())

        job_file.write(job_string.format(**{
            "job_id":       calculation.job_id,
            "molecule":     calculation.molecule.to_xyz(calculation.fragments, calculation.cp).replace("\n", "\\n"),
            "method":       calculation.method,
            "basis":        calculation.basis,
            "num_threads":  settings.get("psi4", "num_threads"),
            "memory":       settings.get("psi4", "memory"),
            "format":       "{}"
        }))
