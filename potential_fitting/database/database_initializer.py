# absolute module imports
from potential_fitting.molecule import parse_training_set_file

# local module imports
from .database import Database

def initialize_database(settings_path, training_set_path, method, basis, cp, *tags):
    """
    Adds energies to be calculated to an existing database.

    Args:
        settings_path       - Local path to the ".ini" file with all relevant settings.
        training_set_path      - Local path to the ".xyz" training set file.
        method              - QM method to use to calculate the energy of these configurations.
        basis               - QM basis to use to calculate the energy of these configurations.
        cp                  - Use counterpoise correction for these configurations?
        tags                - Label this calculation with these tags.

    Returns:
        None.
    """
    
    print("Initializing database from xyz file {} directory into database.".format(training_set_path))

    with open(training_set_path, "r") as training_set_file:
        molecules = parse_training_set_file(filename, settings)

    with Database() as database:
        database.add_calculations(molecules, method, basis, cp, *tags)

    print("Initializing of database successful.")