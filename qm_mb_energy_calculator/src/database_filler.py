import sqlite3
import pickle
import itertools
import psi4
import configparser

import calculator

# this file fills the missing energies in a database
def fill_database(database_name):
    # add .db to the database name if it doesn't already end in .db
    if database_name[-3:] != ".db":
        database_name += ".db"

    # create connection
    connection = sqlite3.connect(database_name)

    # create cursor
    cursor = connection.cursor()

    # get a list of all the rows with missing energies
    cursor.execute("SELECT * FROM Energies WHERE E0='None' OR E1='None' OR E2='None' OR E01='None' OR E12='None' OR E02='None' OR E012='None'")
    rows = cursor.fetchall()

    # loop thru all entries in the database and fill in missing energies
    for energy_row in rows:
        # has id of this molecule
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
                energy = calculator.calc_energy(molecule, combination, config)

                # build energy string to perform insert into table
                entry_string = "E"
                for i in combination:
                    entry_string += str(i)
                # put energy in database
                cursor.execute("UPDATE Energies SET {}=? WHERE ID=? AND model=? AND cp=?".format(entry_string), (energy, ID, model, cp))
    
    # commit changes to database
    connection.commit()
    connection.close()

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

config = configparser.ConfigParser(allow_no_value = False)
config.read("settings.ini")
fill_database("testdb.db")
