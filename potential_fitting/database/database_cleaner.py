# local module imports
from .database import Database
from potential_fitting.utils import SettingsReader
from potential_fitting.molecule import parse_training_set_file


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


def delete_calculations(settings_path, database_config_path, configurations_path, method, basis, cp, *tags):

    molecules = parse_training_set_file(configurations_path, SettingsReader(settings_path))

    with Database(database_config_path) as database:

        database.delete_calculations(molecules, method, basis, cp, *tags)
