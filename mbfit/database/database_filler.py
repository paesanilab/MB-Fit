# external package imports
import sys

# absolute module imports
from mbfit import calculator
from mbfit.calculator import Model
from mbfit.exceptions import LibraryCallError
from mbfit.utils import SettingsReader, files, system

# local module imports
from .database import Database


def fill_database(settings_path, database_config_path, client_name, *tags, calculation_count=sys.maxsize, qm_options={}):
    """
    Loops over uncalculated energies in a database and calculates them.

    Results are submitted to the database in batches. If interrupted
    all results since the last batch will be stuck on "running".
    call clean_database() to set them back to pending.

    Args:
        settings_path       - Local path to the file with all relevant settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        client_name         - Name of the client performing these calculations.
        calculation_count   - Maximum number of calculations to perform. Default is unlimited.
        qm_options           - Dictionary of extra arguments to be passed to the QM code doing the calculation.

    Returns:
        None.
    """

    # open the database
    with Database(database_config_path) as database:

        total_pending = database.count_pending_calculations(*tags)
        system.format_print("Beginning calculations. {} total calculations with tags {} pending in database. Calculating {} of them.".format(total_pending, tags, min(calculation_count, total_pending)),
                bold=True, color=system.Color.YELLOW)


        calc = calculator.get_calculator(settings_path)

        counter = 0
        successes = 0
        failures = 0

        calculation_results = []
        
        for molecule, method, basis, cp, use_cp, frag_indices in database.get_all_calculations(client_name, *tags, calculations_to_do=calculation_count):
            
            counter += 1

            try:
                model = Model(method, basis, use_cp)

                # calculate the missing energy
                energy, log_path = calc.calculate_energy(molecule, model, frag_indices, qm_options=qm_options)
                with open(log_path, "r") as log_file:
                    log_text = log_file.read()
                calculation_results.append((molecule, method, basis, cp, use_cp, frag_indices, True, energy, log_text))
                successes += 1
            
            except LibraryCallError as e:
                if e.log_path is not None:
                    with open(e.log_path, "r") as log_file:
                        log_text = log_file.read()
                    if log_text is "":
                        log_text = "<Log file was empty.>"
                    calculation_results.append((molecule, method, basis, cp, use_cp, frag_indices, False, 0, log_text))
                else:
                    log_text = "<Error occurred without producing log file.>"
                    calculation_results.append((molecule, method, basis, cp, use_cp, frag_indices, False, 0, log_text))
                failures += 1


            if len(calculation_results) >= database.get_batch_size():
                database.set_properties(calculation_results)
                calculation_results = []
                # save changes to the database
                database.save()

            if counter % 10 == 0:
                system.format_print("Performed {} calculations so far. {} Successes and {} Failures so far.".format(counter, successes, failures),
                        italics=True)

        database.set_properties(calculation_results)

        system.format_print("Done! Performed {} calculations. {} Successes and {} Failures. {} calculations with tags {} remain pending in database.".format(counter, successes, failures, total_pending - counter, tags),
                bold=True, color=system.Color.GREEN)


def generate_inputs_from_database(settings_path, database_path):
    """
    Loops over all the uncalculated energies in a database and generates the inputs.
    FIXME Might fail if psi4 is used.

    Can be interrupted, however the energy running when the interrupt occured will be stuck set to running. Once the call is completed, all calculations will be set to running. before running the calculations routine, one should run clean_database().

    Args:
        settings_path       - Local path to the file with all relevant settings information.
        database_path       - Local path to the database file. ".db" will be appending if it does not already end in
                ".db".

    Returns:
        None.
    """
    # open the database
    with Database(database_path) as database:

        print("Starting input creation of database {}".format(database_path))
        # parse settings file
        settings = SettingsReader(settings_path)

        counter = 0

        for calculation in database.missing_energies():

            counter += 1
            print_progress(counter)

            # create input
            calculator.generate_input(calculation.molecule, calculation.fragments, calculation.method + "/" + calculation.basis, calculation.cp, settings)


        print("\nInput creation of database {} successful".format(database_path))

def run_missing_calculations(settings_path, database_path):
    """
    Loops over all the uncalculated energies in a database and runs the calculations.
    FIXME Might fail if psi4 is used. Needs the input generated by generate_inputs_from_database

    Can be interrupted, however the energy running when the interrupt occured will be stuck set to running. Once the call is completed, all calculations will be set to running. before running the routine that retrieves the energy from the outputs, one should run clean_database().

    Args:
        settings_path       - Local path to the file with all relevant settings information.
        database_path       - Local path to the database file. ".db" will be appending if it does not already end in
                ".db".

    Returns:
        None.
    """
    # open the database
    with Database(database_path) as database:

        print("Running missing energies in database {}".format(database_path))
        # parse settings file
        settings = SettingsReader(settings_path)

        counter = 0

        for calculation in database.missing_energies():
            counter += 1
            print_progress(counter)

            # Run this calculation
            calculator.run_calculation(calculation.molecule, calculation.fragments, calculation.method + "/" + calculation.basis, calculation.cp, settings)

def retrieve_energies(settings_path, database_path):
    """
    Loops over all the uncalculated energies in a database and retrieves the energy from the outputs. An error will be raised if there is no output for a given file.
    FIXME Might fail if psi4 is used. 

    Can be interrupted, however the energy running when the interrupt occured will be stuck set to running. If restarting, one should run clean_database().

    Args:
        settings_path       - Local path to the file with all relevant settings information.
        database_path       - Local path to the database file. ".db" will be appending if it does not already end in
                ".db".

    Returns:
        None.
    """
    # open the database
    with Database(database_path) as database:

        print("Retrieving energies to database {}".format(database_path))
        # parse settings file
        settings = SettingsReader(settings_path)

        counter = 0

        for calculation in database.missing_energies():
            counter += 1
            print_progress(counter)

            try:
                # calculate the missing energy
                energy = calculator.retrieve_energy(calculation.molecule, calculation.fragments, calculation.method + "/" + calculation.basis, calculation.cp, settings)
                # update the energy in the database
                database.set_energy(calculation.job_id, energy, "some/log/path")
            except LibraryCallError:
                database.set_failed(calculation.job_id, "some/log/path")

            # save changes to the database
            database.save()

def print_progress(counter):
    """
    Prints the given number to the console. Followed by a newline if the number is divisible by 10.

    Args:
        counter             - The number to print.

    Returns:
        None.
    """

    if counter % 10 == 0:
        s = "Beginning calculation number {:6d}.".format(counter)
        s += "\n"
        print(s, end="", flush=True)
