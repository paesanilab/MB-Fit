# local module imports
from .database import Database

def clean_database(settings_path):
    """
    Sets all running calculations back to pending in the given database.

    Args:
        settings_path       - Local path to ".ini" file containing all relevent settings.

    Returns:
        None.
    """

    with Database() as database:

        database.clean()
