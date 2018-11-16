# absolute module imports
from potential_fitting.molecule import Molecule

# local module imports
from .database import Database, Calculation

def read_job(database_path, job_path, job_log_path):
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

    with open(job_path, "r") as job_file:

        # parse the job id
        job_id = int(job_file.readline().split()[1])
        
        output_line = job_file.readline()[:-1]

        if output_line == "Failure":
            success = False
        else:
            # parse the energy2
            energy = float(output_line.split()[1])
            success = True
    
    # open the database
    with Database(database_path) as database:

        if success:
            database.set_energy(job_id, energy, job_log_path)
        else:
            database.set_failed(job_id, "failed", job_log_path)
