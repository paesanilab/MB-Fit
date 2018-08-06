import sys
import configparser
from database import Database
def make_all_jobs(settings, database_name, job_dir):
    """
    Makes all the jobs that still need to be performed in this database

    Args:
        settings - .ini file with relevent settings
        database_name - file path to where the database is stored
        job_dir - directory to place the jobs in

    Returns:
        None
    """
    
    # add .db to the database name if it doesn't already end in .db
    if database_name[-3:] != ".db":
        database_name += ".db"

    # open the database
    database = Database(database_name)

    for calculation in database.missing_energies():
        write_job(settings, calculation, job_dir)

    # commit changes to database
    database.save()
    # close the database
    database.close()

def make_job(settings, database_name, job_dir):   
    """
    Makes a single job that still needs to be performed in this database

    Args:
        settings - .ini file with relevent settings
        database_name - file path to where the database is stored
        job_dir - directory to place the job in

    Returns:
        None
    """
    
    # add .db to the database name if it doesn't already end in .db
    if database_name[-3:] != ".db":
        database_name += ".db"

    # open the database
    database = Database(database_name)

    write_job(settings, database.get_missing_energy(), job_dir)

    # commit changes to database
    database.save()
    # close the database
    database.close()

def write_job(settings, calculation, job_dir):
    """
    Makes a job from a specific Calculation

    file will be located at job_dir/job_<id>.py
    
    Args:
        settings - .ini file with relevent settings
        calculation - the Calculation object (see database.py) with the information needed to make a job
        job_dir - directory to place the job in

    Returns:
        None
    
    """
    # parse settings file
    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)

    with open(job_dir + "/job_{}.py".format(calculation.job_id), "w") as job_file, open("job_template.py", "r") as job_template:
        job_string = "".join(job_template.readlines())
        
        job_file.write(job_string.format(**{
            "job_id":       calculation.job_id,
            "molecule":     calculation.molecule.to_xyz(calculation.fragments, calculation.cp).replace("\n", "\\n"),
            "method":       calculation.method,
            "basis":        calculation.basis,
            "num_threads":  config["psi4"]["num_threads"],
            "memory":       config["psi4"]["memory"],
            "format":       "{}"
        }))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python database_job_maker.py <settings> <database_name> <job_dir>")
        sys.exit(1)
    make_all_jobs(sys.argv[1], sys.argv[2], sys.argv[3])
