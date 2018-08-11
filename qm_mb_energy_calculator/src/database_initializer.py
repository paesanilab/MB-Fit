import os, sys

import sqlite3, configparser
from molecule_parser import xyz_to_molecules
from database import Database

"""
initializes a database from the config files in a directory
"""
def initialize_database(settings, database_name, directory):
    """
    Initializes a new database or adds energies to be calculated to an existing one.

    Args:
        settings - the .ini file used to initialize the database
        database_name - the name of the file to save the database in
        directory - the directory of the config.xyz files to add to the database

    Returns:
        None
    """

    # add .db to database name if it doesn't already end in .db 
    if database_name[-3:] != ".db":
        print("Database name \"{}\" does not end in database suffix \".db\". Automatically adding \".db\" to end of database name.".format(database_name))
        database_name += ".db"

    # Create the database
    database = Database(database_name)
    database.create()
    
    print("Initializing database from xyz files in {} directory into database {}".format(directory, database_name))
    # get a list of all the files in the directory, or just the filename if directory is just a single file.
    filenames = get_filenames(directory) if os.path.isdir(directory) else [directory]

    # parse settings.ini file
    config = configparser.ConfigParser(allow_no_value=False)

    # See if a config file already exists; if not, generate one
    config.read(settings)
    
    # Read values from the config file
    method = config["energy_calculator"]["method"]
    basis = config["energy_calculator"]["basis"]
    cp = config["energy_calculator"].getboolean("cp")
    tag = config["molecule"]["tag"]

    # loop thru all files in directory
    for filename in filenames:

        # if the filename does not end in .xyz, then skip it
        if filename[-4:] != ".xyz":
            continue

        # open the file
        xyz_file = open(filename, "r")

        # get list of all molecules in file
        molecules = xyz_to_molecules(xyz_file, config)
        
        # if this file contains omptimized geometries, tell the database so
        if filename[-8:] == ".opt.xyz":

            # for each molecule in the file
            for molecule in molecules:

                # add this molecule to the database, optimized flag set to true
                database.add_calculation(molecule, method, basis, cp, tag, True)

        else: 

            # for each molecule in the file
            for molecule in molecules:

                # add this molecule to the database, optimized flag set to false
                database.add_calculation(molecule, method, basis, cp, tag, False)

    # save and close the database
    database.save()
    database.close()
    
    print("Initializing of database {} successful".format(database_name))


def get_filenames(directory):
    """
    gets list of all filenames in a directory

    Args:
        directory - directory to get filenames from

    Returns:
        list of all filenames in that directory, and subdirectories
    """
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
