import sys, os
import sqlite3
import calculator
from database import Database

from potential_fitting.exceptions import LibraryCallError
from potential_fitting.utils import SettingsReader

def fill_database(settings_file, database_name):
    """
    Calculates all the pending energies in a database

    Args:
        settings_file   - the file with all relevant settings information
        database_name   - the database file

    Returns:
        None
    """

    # open the database
    with Database(database_name) as database:

        print("Filling database {}".format(database_name))
        # parse settings file
        settings = SettingsReader(settings_file)

        counter = 0
        
        for calculation in database.missing_energies():
            
            counter += 1
            print_progress(counter)

            try:
                # calculate the missing energy
                energy = calculator.calculate_energy(calculation.molecule, calculation.fragments, calculation.method + "/" + calculation.basis, calculation.cp, settings)
                # update the energy in the database
                database.set_energy(calculation.job_id, energy, "some/log/path")
            except LibraryCallError:
                database.set_failed(calculation.job_id, "failed", "some/log/path")
            
            # save changes to the database
            database.save()

        print("\nFilling of database {} successful".format(database_name))

def print_progress(counter):
    s = "{:6d}".format(counter)
    if counter % 10 == 0:
       s += "\n" 
    print(s, end="", flush=True)
