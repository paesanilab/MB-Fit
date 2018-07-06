import sys
import sqlite3
import pickle
from database import Database

def generate_fitting_input(settings, database_name, output_path):
    # add .db to database name if it doesn't already end in .db 
    if database_name[-3:] != ".db":
        print("Database name \"{}\" does not end in database suffix \".db\". Automatically adding \".db\" to end of database name.".format(database_name))
        database_name += ".db"
    
    database = Database(database_name)

    print("Creating a fitting input file from database {} into file {}".format(database_name, output_path))

    molecule_energy_pairs = database.get_complete_energies()

    # open file for writing
    output = open(output_path, "w");

    for molecule_energy_pair in molecule_energy_pairs:
        molecule = molecule_energy_pair[0]
        energies = molecule_energy_pair[1]
        output.write(str(molecule.get_num_atoms()) + "\n")
        
        for energy in energies:
            output.write(str(energy) + " ")

        output.write("\n")

        output.write(molecule.to_xyz() + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Incorrect number of arguments");
        print("Usage: python database_reader.py <settings_file> <database_name> <output_path>")
        sys.exit(1)   
    generate_fitting_input(sys.argv[1], sys.argv[2], sys.argv[3])
