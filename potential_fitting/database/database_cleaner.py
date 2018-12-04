# local module imports
from .database import Database

def clean_database(settings_path, database_path):
    """
    Sets all running calculations back to pending in the given database.

    Args:
        settings_path       - Local path to ".ini" file containing all relevent settings.
        database_path       - Local path to the file with the database. ".db" will be appended if it does not already
                end in ".db".

    Returns:
        None.
    """

    with Database(database_path) as database:

        database.clean()
