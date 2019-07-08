# local module imports
from .database import Database

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