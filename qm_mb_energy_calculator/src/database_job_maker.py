import sys
import configparser
from database import Database
def make_all_jobs(settings, database_name):
    
    # add .db to the database name if it doesn't already end in .db
    if database_name[-3:] != ".db":
        database_name += ".db"

    # open the database
    database = Database(database_name)

    for calculation in database.missing_energies():
        write_job(settings, calculation)

    # commit changes to database
    database.save()
    # close the database
    database.close()

def make_job(settings, database_name):   
    # add .db to the database name if it doesn't already end in .db
    if database_name[-3:] != ".db":
        database_name += ".db"

    # open the database
    database = Database(database_name)

    write_job(settings, database.get_missing_energy())

    # commit changes to database
    database.save()
    # close the database
    database.close()

def write_job(settings, calculation):
    # parse settings file
    config = configparser.SafeConfigParser(allow_no_value=False)
    config.read(settings)

    # header for job
    print("Job: {}".format(calculation.job_id))
    
    # print the molecule's xyz
    print(calculation.molecule.to_xyz(calculation.fragments, calculation.cp))

    # print the method, basis, and cp
    print("Method: {}".format(calculation.method))
    print("Basis: {}".format(calculation.basis))

    # other settings stuff
    print("Code: {}".format(config["energy_calculator"]["code"]))
    print("Memory: {}".format(config["psi4"]["memory"]))
    print("Threads: {}".format(config["psi4"]["num_threads"]))

    print()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python database_job_maker.py <settings> <database_name>")
        sys.exit(1)
    make_all_jobs(sys.argv[1], sys.argv[2])
