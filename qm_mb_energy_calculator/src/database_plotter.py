import sys
import sqlite3
from database import Database

def plot_difference(database_name, method1, method2, basis1, basis2, cp1, cp2, energy):
    """
    Takes in two different sets of method/basis/cp and prints the average difference in energy for the same molecules between them
    """
    
    # add .db to the database name if it doesn't already end in .db
    if database_name[-3:] != ".db":
        print("Database name \"{}\" does not end in database suffix \".db\". Automatically adding \".db\" to end of database name.".format(database_name))
        database_name += ".db"

    database = Database(database_name)
    
    # get array of pairs of energies computed in the two different ways
    energy_pairs = database.get_comparison_energies(method1, method2, basis1, basis2, cp1, cp2, energy)

    differences = []    
    for energy_pair in energy_pairs:
        # make sure both energies are numbers and not "N/A" MIGHT BE UNNEEDED?
        if not energy_pair[0] == "N/A" and not energy_pair[1] == "N/A":
            differences.append(energy_pair[0] - energy_pair[1])
    try:
        print("Average Difference in Energy: {:3.3e}".format(sum(differences) / len(differences)))
    except ZeroDivisionError as error:
        raise ValueError("No Matching energies for both categories for those methods and basises.") from error

if __name__ == "__main__":
    if sys.argv[1] == "difference":
        plot_difference(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9])
