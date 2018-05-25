import sqlite3
import os
from molecule_parser import xyz_to_molecules
import database
import pickle

"""
initializes a database from the config files in a directory
"""
def initialize_database(database_name, directory, tag, model, cp):
    # add .db to database name if it doesn't already end in .db 
    if database_name[-3:] != ".db":
        database_name += ".db"

    # create connection
    connection = sqlite3.connect(database_name)

    # create cursor
    cursor = connection.cursor()

    # create the tables
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Configs(ID text, config BLOB, natoms INT,
        nfrags INT, tag TEXT)
        """)
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Energies(ID TEXT, model TEXT, cp INT,
        E0 REAL, E1 REAL, E2 REAL, E01 REAL, E02 REAL, E12 REAL, E012 REAL, Enb BLOB)
        """)
    
    if not os.path.isdir(directory):
        raise Exception("Directory {} does not exist".format(directory))
    filenames = get_filenames(directory)

    # loop thru all files in directory
    for filename in filenames:
        # open the file
        f = open(filename, "r")
        # get list of all molecules in file
        molecules = xyz_to_molecules(f)
        # for each molecule in the file
        for molecule in molecules:
            hash_id = molecule.get_SHA1()
            # check if molecule already has rows in the table
            cursor.execute("select ID from Configs where ID=?", (hash_id,))
            config_row = cursor.fetchone()
            cursor.execute("select ID from Energies where ID=?", (hash_id,))
            energies_row = cursor.fetchone()
            
            # create a new row in the Configs table if such a row did not already exist
            if config_row is None:
                pick = pickle.dumps(molecule)
                config = pick.hex() # turn byte array into hex string
                # convert hex string into byte array bytes.fromhex(config)
                cursor.execute("INSERT INTO Configs (ID, config, natoms, nfrags, tag) VALUES ('{}', '{}', '{}', '{}', '{}')".format(hash_id, config, molecule.get_num_atoms(), molecule.get_num_fragments(), tag))

            # create a new row in the Energies table if such a row did not already exist
            if energies_row is None:
                cursor.execute("INSERT INTO Energies (ID, model, cp, E0, E1, E2, E01, E02, E12, E012, Enb) VALUES ('{}', '{}', '{}', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None')".format(hash_id, model, cp))
    
    connection.commit()
    connection.close()


"""
gets list of all filenames in a director
"""
def get_filenames(directory):
    filenames = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            filenames.append(root + "/" + filename);
    return filenames;
    

initialize_database("testdb.db", "config_files", "someTag", "someModel", "cpT/F")
