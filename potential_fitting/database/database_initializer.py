# absolute module imports
from potential_fitting.molecule import parse_training_set_file
from potential_fitting.utils import SettingsReader

# local module imports
from .database import Database

def initialize_database(settings_path, database_config_path, training_set_path, method, basis, cp, *tags, optimized = False):
    """
    Adds energies to be calculated to an existing database.

    Args:
        settings_path       - Local path to the ".ini" file with all relevant settings.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        training_set_path      - Local path to the ".xyz" training set file.
        method              - QM method to use to calculate the energy of these configurations.
        basis               - QM basis to use to calculate the energy of these configurations.
        cp                  - Use counterpoise correction for these configurations?
        tags                - Label this calculation with these tags.
        optimized 			- Are these configurations optimized geometries? Default is False.

    Returns:
        None.
    """
    
    print("Initializing database from xyz file {} directory into database.".format(training_set_path))

    molecules = parse_training_set_file(training_set_path, SettingsReader(settings_path))

    with Database(database_config_path) as database:
        database.add_calculations(molecules, method, basis, cp, *tags, optimized = optimized)

    print("Initializing of database successful.")