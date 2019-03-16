# absolute module imports
from potential_fitting.molecule import xyz_to_molecules

# local module imports
from .database import Database

def initialize_database(settings_path, training_set_path, method, basis, cp, *tags):
    """
    Initializes a new database or adds energies to be calculated to an existing one.

    Args:
        settings_path       - Local path to the ".ini" file with all relevant settings.
        database_path       - Local path to the database file. ".db" will be appended if it does not already end in
                ".db".
        directory_path      - Local path to the directory of ".xyz" files or a single ".xyz" file to add to the
                database. If a directory is specified, all non ".xyz" or ".xyz.opt" files will be ignored"
        tags                - Label this calculation with these tags.

    Returns:
        None.
    """
    
    print("Initializing database from xyz file {} directory into database.".format(training_set_path))

    with open(training_set_path, "r") as training_set_file:
        molecules = xyz_to_molecules(filename, settings)

    print("Training set size: {}".format(len(molecules)))

    with Database() as database:
        database.add_calculations(molecules, method, basis, cp, *tags)

    print("Initializing of database successful.")