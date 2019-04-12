# local module imports
from .database import Database

def clean_database(database_config_path, settings_path):
    """
    Sets all running calculations back to pending in the given database.

    Args:
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        settings_path       - Local path to ".ini" file containing all relevent settings.

    Returns:
        None.
    """

    with Database(config_path) as database:

        database.clean()
