import sys
from molecule import Molecule
from database import Database, Calculation

def read_job(database_name, job_path, job_log_path):
    with open(job_path, "r") as job_file:

        # skip the "Job" line
        job_id = int(job_file.readline().split()[1])

        # parse the energy
        energy = float(job_file.readlines()[-1].split()[1])

    # add .db to the database name if it doesn't already end in .db
    if database_name[-3:] != ".db":
        database_name += ".db"

    # open the database
    database = Database(database_name)

    database.set_energy(job_id, energy, job_log_path)

    # commit changes to database
    database.save()
    # close the database
    database.close()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Incorrect number of arguments");
        print("Usage: python database_filler.py <database_name> <job_path> <job_log_path>")
    read_job(sys.argv[1], sys.argv[2], sys.argv[3])
