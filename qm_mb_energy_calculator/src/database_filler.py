import sqlite3
import pickle
import itertools
import psi4

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
        cursor.execute("SELECT * from Configs WHERE ID=?", (ID,))
        config_row = cursor.fetchone()

        # load molecule object
        molecule = pickle.loads(bytes.fromhex(config_row[1]))

        # loop thru all the energies in the row
        for index, energy_entry in enumerate(energy_row[3:10]):
            # check if energy entry is missing
            if energy_entry == "None":
                combination = get_combination_from_index(index)
                psi4_string = molecule.to_xyz(combination)
                psi4_mol = psi4.core.Molecule.create_molecule_from_string(psi4_string)
                psi4_mol.update_geometry()
                psi4.set_num_threads(1)

                # calculate energy using psi4
                energy = psi4.energy("HF/STO-3G", molecule=psi4_mol)

                # build energy string to perform insert into table
                entry_string = "E"
                for i in combination:
                    entry_string += str(i)
                print("TEST: {}".format(energy)); 
                # put energy in database
                cursor.execute("UPDATE Energies SET {}=? WHERE ID=? AND model=? AND cp=?".format(entry_string), (energy, ID, model, cp))
    
    connection.commit();
    # temporary tests
    cursor.execute("SELECT * from Energies")
    print(cursor.fetchall())
    connection.close();


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
fill_database("testdb.db")
