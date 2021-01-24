# local module imports
from .database import Database
from mbfit.utils import SettingsReader
from mbfit.molecule import parse_training_set_file


def clean_database(settings_path, database_config_path, *tags):
    """
    Sets all dispatched calculations back to pending in the given database.

    Args:
        settings_path       - Local path to ".ini" file containing all relevent settings.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        tags                - Reset calculations with one of these tags.

    Returns:
        None.
    """

    with Database(database_config_path) as database:

        database.reset_dispatched(*tags)


def reset_database(settings_path, database_config_path, *tags):
    """
    Sets all failed calculations back to pending in the given database.

    Args:
        settings_path       - Local path to ".ini" file containing all relevent settings.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        tags                - Reset calculations with one of these tags.

    Returns:
        None.
    """

    with Database(database_config_path) as database:

        database.reset_failed(*tags)


def delete_calculations(settings_path, database_config_path, configurations_path, method, basis, cp, *tags, delete_complete_calculations = False):
    """
    Removes the specified tags from any calculations in the database that matches one of the molecules in
    the configurations file and the given method, basis, and cp.

    Will never delete completed or dispatched energies unless delete_complete_calculations is True, only remove tags from them.

    Will fully delete uncomplete calculations from the database.

    Args:
        settings_path       - Local path to ".ini" file containing all relevent settings.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        configurations_path - '.xyz' file. Remove tags from calculations involving these molecules.
        method  - Remove tags from calculations with this method.
        basis   - Remove tags from calculations with this basis.
        cp      - Remove tags from calculations with this cp.
        tags    - The tags to remove.
        delete_complete_calculations - If True, delete calculations even if their energy is already
                calculated.

    Returns:
        None.
    """
    molecules = parse_training_set_file(configurations_path, SettingsReader(settings_path))

    with Database(database_config_path) as database:

        database.delete_calculations(molecules, method, basis, cp, *tags, delete_complete_calculations=delete_complete_calculations)


def delete_all_calculations(settings_path, database_config_path, molecule_name, method, basis, cp, *tags, delete_complete_calculations = False):
    """
    Removes the specified tags from any calculations in the database that matches one of the molecules in
    the configurations file and the given method, basis, and cp.

    Will never delete completed or dispatched energies unless delete_complete_calculations is True, only remove tags from them.

    Will fully delete uncomplete calculations from the database.

    Args:
        settings_path       - Local path to ".ini" file containing all relevent settings.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        molecule_name - Remove tags from calculations involving molecules of this type.
        method  - Remove tags from calculations with this method.
        basis   - Remove tags from calculations with this basis.
        cp      - Remove tags from calculations with this cp.
        tags    - The tags to remove.
        delete_complete_calculations - If True, delete calculations even if their energy is already
                calculated.

    Returns:
        None.
    """

    with Database(database_config_path) as database:

        database.delete_all_calculations(molecule_name, method, basis, cp, *tags, delete_complete_calculations=delete_complete_calculations)

