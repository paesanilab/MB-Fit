# external package imports
import sys

# absolute module imports
from potential_fitting import calculator
from potential_fitting.exceptions import LibraryCallError
from potential_fitting.utils import SettingsReader, files

# local module imports
from .database import Database


def fill_database(settings_path, client_name, calculation_count=sys.maxsize):
    """
    Loops over uncalculated energies in a database and calculates them.

    Results are submitted to the database in batches. If interrupted
    all results since the last batch will be stuck on "running".
    call clean_database() to set them back to pending.

    Args:
        settings_path       - Local path to the file with all relevant settings information.
        client_name         - Name of the client performing these calculations.
        calculation_count   - Maximum number of calculations to perform. Default is unlimited.

    Returns:
        None.
    """

    # open the database
    with Database() as database:

        print("Filling database.")

        # parse settings file
        settings = SettingsReader(settings_path)

        counter = 0

        calculation_results = []
        
        for molecule, method, basis, cp, use_cp, frag_indices in database.get_all_calculations(client_name, calculations_to_do=calculation_count):
            
            counter += 1
            print_progress(counter)

            try:

                # calculate the missing energy
                energy = calculator.calculate_energy(molecule, frag_indices, method + "/" + basis, use_cp, settings)
                calculation_results.append((molecule, method, basis, cp, use_cp, frag_indices, True, energy, "log_text"))
            
            except LibraryCallError:

                calculation_results.append((molecule, method, basis, cp, use_cp, frag_indices, False, 0, "log_text"))

            if len(calculation_results) >= database.get_batch_size():
                database.set_properties(calculation_results)
                # save changes to the database
                database.save()

        database.set_properties(calculation_results)

        print("\nFilling of database successful")


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
