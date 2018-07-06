import sqlite3
import os
from molecule_parser import xyz_to_molecules
import pickle
import configparser
import sys
from database import Database

"""
initializes a database from the config files in a directory
"""
def initialize_database(settings, database_name, directory):
    # add .db to database name if it doesn't already end in .db 
    if database_name[-3:] != ".db":
        print("Database name \"{}\" does not end in database suffix \".db\". Automatically adding \".db\" to end of database name.".format(database_name))
        database_name += ".db"

    database = Database(database_name)

    if not os.path.isdir(directory):
        raise ValueError("{} is not a directory. \n Terminating database initialization.".format(directory))
        sys.exit(1)
    
    database.create()
    
    # list of all filenames to parse
    filenames = []

    print("Initializing database from xyz files in {} directory into database {}".format(directory, database_name))
    filenames = get_filenames(directory)

    # parse settings.ini file
    config = configparser.SafeConfigParser(allow_no_value=False)

    # See if a config file already exists; if not, generate one
    config.read(settings)
    
    # Read values from the config file
    method = config["energy_calculator"]["method"]
    basis = config["energy_calculator"]["basis"]
    cp = config["energy_calculator"]["cp"]
    tag = config["molecule"]["tag"]

    # loop thru all files in directory
    for filename in filenames:
        if filename[-4:] != ".xyz":
            continue
        # open the file
        f = open(filename, "r")
        # get list of all molecules in file
        molecules = xyz_to_molecules(f, config)
        # for each molecule in the file
        for molecule in molecules:
            # add this molecule to the database
            database.add_molecule(molecule, method, basis, cp, tag)

    database.save()
    database.close()
    
    print("Initializing of database {} successful".format(database_name))


"""
gets list of all filenames in a directory
"""
def get_filenames(directory):
    filenames = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            filenames.append(root + "/" + filename);
    return filenames;

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Incorrect number of arguments");
        print("Usage: python database_initializer.py <settings_file> <database_name> <config_directory>")
        sys.exit(1)   
    initialize_database(sys.argv[1], sys.argv[2], sys.argv[3])
