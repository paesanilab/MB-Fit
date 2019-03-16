# external package imports
import sys, os, sqlite3

# absolute module imports
from potential_fitting import calculator
from potential_fitting.exceptions import LibraryCallError
from potential_fitting.utils import SettingsReader

# local module imports
from .database import Database

def fill_database(settings_path, client_name):
    """
    Loops over all the uncalculated energies in a database and calculates them.

    Can be interrupted, however the energy running when the interrupt occured will be stuck set to running. Call
    clean_database() to reset it to pending.

    Args:
        settings_path       - Local path to the file with all relevant settings information.
        database_path       - Local path to the database file. ".db" will be appending if it does not already end in
                ".db".

    Returns:
        None.
    """

    # open the database
    with Database() as database:

        print("Filling database.")

        # parse settings file
        settings = SettingsReader(settings_path)

        counter = 0
        
        for molecule, method, basis, cp, frag_indices in database.get_all_calculations(client_name):
            
            counter += 1
            print_progress(counter)

            try:
                # calculate the missing energy
                energy = calculator.calculate_energy(molecule, frag_indices, method + "/" + basis, cp, settings)
                # update the energy in the database
                database.set_energy(calculation.job_id, energy, "some/log/path")
            except LibraryCallError:
                database.set_failed(calculation.job_id, "some/log/path")
            
            # save changes to the database
            database.save()

        print("\nFilling of database {} successful".format(database_path))

def print_progress(counter):
    """
    Prints the given number to the console. Followed by a newline if the number is divisible by 10.

    Args:
        counter             - The number to print.
    
    Returns:
        None.
    """

    s = "{:6d}".format(counter)
    if counter % 10 == 0:
       s += "\n" 
    print(s, end="", flush=True)
