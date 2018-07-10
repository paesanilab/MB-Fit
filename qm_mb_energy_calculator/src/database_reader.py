import sys
import os

# this is messy, I hope there is a better way to do this!
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import sqlite3
import pickle
from database import Database
import constants

def generate_fitting_input(settings, database_name, output_path):
    """
    Creates a training set file from the calculated energies in a database
    
    Currently ONLY WORKS FOR MONOMERS
    """
    
    # add .db to database name if it doesn't already end in .db 
    if database_name[-3:] != ".db":
        print("Database name \"{}\" does not end in database suffix \".db\". Automatically adding \".db\" to end of database name.".format(database_name))
        database_name += ".db"
    
    database = Database(database_name)

    print("Creating a fitting input file from database {} into file {}".format(database_name, output_path))

    # get list of all [molecule, energies] pairs calculated in the database where energies is a list: [E0, E1, E2, E01, E12, E02, E012] with N/A energies missing
    molecule_energy_pairs = database.get_complete_energies()

    # if there are no calculated energies, error and exit
    if len(molecule_energy_pairs) == 0:
        print("No completed energy entries to generate a training set.")
        database.close()
        exit(1)
    
    # temporary code only works for 1b
    
    # find the lowest energy. In future, lowest energy should be energy of optimized geometry, but for now, it is assumed the database contains the optimized geometry (or something close enough)
    min_energy = molecule_energy_pairs[0][1][0]
    for molecule_energy_pair in molecule_energy_pairs:
        min_energy = min_energy if min_energy < molecule_energy_pair[1][0] else molecule_energy_pair[1][0]

    # open output file for writing
    output = open(output_path, "w")

    for molecule_energy_pair in molecule_energy_pairs:
        molecule = molecule_energy_pair[0]
        energies = molecule_energy_pair[1]

        # write the number of atoms to the output file
        output.write(str(molecule.get_num_atoms()) + "\n")
        
        # write the energies to the output file
        for energy in energies:
            output.write(str((float(energy) - float(min_energy)) * constants.au_to_kcal) + " ") # covert Hartrees to kcal/mol

        output.write("\n")

        # write the molecule's atoms' coordinates to the xyz file
        output.write(molecule.to_xyz() + "\n")

    database.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Incorrect number of arguments");
        print("Usage: python database_reader.py <settings_file> <database_name> <output_path>")
        sys.exit(1)   
    generate_fitting_input(sys.argv[1], sys.argv[2], sys.argv[3])
