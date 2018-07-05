import sys
import sqlite3
import pickle

def generate_fitting_input(settings, database_name, output_path):
    # add .db to database name if it doesn't already end in .db 
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
        print("{} exists but is not a valid database file. \n Terminating database initialization.".format(database_name))
        sys.exit(1)

    print("Creating a fitting input file from database {} into file {}".format(database_name, output_path))

    # get a list of all the rows with filled energies
    cursor.execute("SELECT * FROM Energies WHERE NOT (E0='None' OR E1='None' OR E2='None' OR E01='None' OR E12='None' OR E02='None' OR E012='None')")
    rows = cursor.fetchall()

    # open file for writing
    output = open(output_path, "w");
    
    # loop thru all the rows with filled energies
    for energy_row in rows:
        # add each row to output file in proper format
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

        num_frags = molecule.get_num_fragments()
        num_atoms = molecule.get_num_atoms()

        if num_frags == 3:
            # write number of atoms
            output.write(str(num_atoms) + "\n")
            # write energies in comment line
            for i in range(3, 10):
                output.write(str(energy_row[i]) + " ")
            output.write("\n")
            # write xyz coordinates
            output.write(molecule.to_xyz() + "\n")
            pass
        elif num_frags == 2:
            # write number of atoms
            output.write(str(num_atoms) + "\n")
            # write energies in comment line
            for i in [3, 4, 6]:
                output.write(str(energy_row[i]) + " ")
            output.write("\n")
            # write xyz coordinates
            output.write(molecule.to_xyz() + "\n")
            pass
        elif num_frags == 1:
            # write number of atoms
            output.write(str(num_atoms) + "\n")
            # write energies in comment line
            output.write(str(energy_row[3]) + " ")
            output.write("\n")
            # write xyz coordinates
            output.write(molecule.to_xyz() + "\n")
            pass
        else:
            print("Unsupported Number of fragments {}. Supported values are 1,2, and 3.".format(fragment_count))
            

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Incorrect number of arguments");
        print("Usage: python database_reader.py <settings_file> <database_name> <output_path>")
        sys.exit(1)   
    generate_fitting_input(sys.argv[1], sys.argv[2], sys.argv[3])
