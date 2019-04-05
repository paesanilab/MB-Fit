# absolute module imports
from potential_fitting.molecule import Molecule, os
from glob import glob

# local module imports
from .database import Database

def read_all_jobs(job_dir):
    calculation_results = []
    for directory in glob("job_[0-9]+"):
        calculation_results.append(read_job(directory + "/output.dat", directory + "/output.log"))

        i = 1

        job_dir = "job_{format}_done".format(i)

        while os.path.isfile(job_dir):
            i += 1

            job_dir = "job_{format}_done".format(i)

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

    with open(job_dat_path, "r") as job_file:

        molecule =

        if output_line == "Failure":
            success = False
        else:
            # parse the energy2
            energy = float(output_line.split()[1])
            success = True


    return molecule, method, basis, cp, frag_indices, True, energy, "log_text"

