import sys, os
import sqlite3
import calculator
from database import Database

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
from exceptions import LibraryCallError
import settings_reader

def fill_database(settings_file, database_name, directory = "unused"): # argument is unused, but I haven't removed yet because it will break a lot of code
    """
    Walks through a database and calculates all missing energies
    
    settings_file is the .ini file to use to fill this database

    database_name is file the database is stored in

    directory is unued. # should probably be removed.

    """
    # add .db to the database name if it doesn't already end in .db
    if database_name[-3:] != ".db":
        print("Database name \"{}\" does not end in database suffix \".db\". Automatically adding \".db\" to end of database name.".format(database_name))
        database_name += ".db"

    with Database(database_name) as database:

        print("Filling database {}".format(database_name))

        # parse settings file
        settings = settings_reader.SettingsReader(settings_file)
        
        for calculation in database.missing_energies():

            try:
                # calculate the missing energy
                energy = calculator.calculate_energy(calculation.molecule, calculation.fragments, calculation.method + "/" + calculation.basis, calculation.cp, settings)
                # update the energy in the database
                database.set_energy(calculation.job_id, energy, "some/log/path")
            except LibraryCallError:
                database.set_failed(calculation.job_id, "failed", "some/log/path")
            
            # save changes to the database
            database.save()

        print("Filling of database {} successful".format(database_name))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Incorrect number of arguments");
        print("Usage: python database_filler.py <settings_file> <database_name> <directory>")
        exit(1)
    fill_database(sys.argv[1], sys.argv[2], sys.argv[3])
