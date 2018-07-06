import sys
import sqlite3

def plot_difference(database_name, model1, model2, cp1, cp2, energy):
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
        raise ValueError("{} exists but is not a valid database file. \n Terminating database initialization.".format(database_name))
        sys.exit(1)

    # get rows pertaining to category 2
    cursor.execute("SELECT ID FROM Energies WHERE model=? AND cp=?", (model1, cp1))
    category1_IDs = cursor.fetchall()
    cursor.execute("SELECT {} FROM Energies WHERE model=? AND cp=?".format(energy), (model1, cp1))
    category1_energies = cursor.fetchall()

    # get rows pertaining to category 2
    cursor.execute("SELECT ID FROM Energies WHERE model=? AND cp=?", (model2, cp2))
    category2_IDs = cursor.fetchall()
    cursor.execute("SELECT {} FROM Energies WHERE model=? AND cp=?".format(energy), (model2, cp2))
    category2_energies = cursor.fetchall()

    # make list of all ids
    IDs = category1_IDs.copy()
    for ID in category2_IDs:
        if ID not in category1_IDs:
            IDs.append(ID)

    # graph the difference
    differences = []
    for ID in IDs:
        try:
            category1_energy = category1_energies[category1_IDs.index(ID)][0]
            category2_energy = category2_energies[category2_IDs.index(ID)][0]
            
            difference = category1_energy - category2_energy

            differences.append(difference)

            print("Difference in Energy: {:3.3e}".format(difference))
        except ValueError:
            pass
    try:
        print("Average Difference in Energy: {:3.3e}".format(sum(differences) / len(differences)))
    except ZeroDivisionError as error:
        raise ValueError("No Matching energies for both categories for those methods and basises.") from error

if __name__ == "__main__":
    if sys.argv[1] == "difference":
        plot_difference(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
