import sqlite3
import pickle
import itertools
import psi4
import configparser
import sys

import calculator

# this file fills the missing energies in a database
def fill_database(settings, database_name, directory):
    # add .db to the database name if it doesn't already end in .db
    if database_name[-3:] != ".db":
        print("Database name \"{}\" does not end in database suffix \".db\". Automatically adding \".db\" to end of database name.".format(database_name))
        database_name += ".db"
    
    # create connection
    connection = sqlite3.connect(database_name)

    # create cursor
    cursor = connection.cursor()

    try:
        cursor.execute("PRAGMA table_info('schema_version')");
    except:
        print("{} exists but is not a valid database file. \n Terminating database filling.".format(database_name))
        sys.exit(1)

    print("Filling database {}".format(database_name))
      

    # get a list of all the rows with missing energies
    cursor.execute("SELECT * FROM Energies WHERE E0='None' OR E1='None' OR E2='None' OR E01='None' OR E12='None' OR E02='None' OR E012='None'")
    rows = cursor.fetchall()

    # parse settings.ini file
    config = configparser.SafeConfigParser(allow_no_value=False)

    # See if a config file already exists; if not, generate one
    config.read(settings)

    # loop thru all entries in the database and fill in missing energies
    for energy_row in rows:
        # hash id of this molecule
        ID = energy_row[0]
        # model to use to calc energy of this molecule
        model = energy_row[1]
        # wheter to use cp on calcs
        cp = energy_row[2]
        # fetch the Configs row
        cursor.execute("SELECT * FROM Configs WHERE ID=?", (ID,))
        config_row = cursor.fetchone()

        # load molecule object
        molecule = pickle.loads(bytes.fromhex(config_row[1]))

        # loop thru all the energies in the row
        for index, energy_entry in enumerate(energy_row[3:10]):
          
            # check if energy entry is missing
            if energy_entry == "None":
                combination = get_combination_from_index(index)

                # Call the calculator function
                energy = calculator.calculate_energy(molecule, combination, model, cp, config)

                # build energy string to perform insert into table
                entry_string = "E"
                for i in combination:
                    entry_string += str(i)
                # put energy in database
                cursor.execute("UPDATE Energies SET {}=? WHERE ID=? AND model=? AND cp=?".format(entry_string), (energy, ID, model, cp))
    
    # commit changes to database
    connection.commit()
    connection.close()

    print("Filling of database {} successful".format(database_name))

# generates a combination from an index
def get_combination_from_index(index):
    return [
            [0],
            [1],
            [2],
            [0, 1],
            [0, 2],
            [1, 2],
            [0, 1, 2]
           ][index]

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Incorrect number of arguments");
        print("Usage: python database_filler.py <settings_file> <database_name> <directory>")
    fill_database(sys.argv[1], sys.argv[2], sys.argv[3])
