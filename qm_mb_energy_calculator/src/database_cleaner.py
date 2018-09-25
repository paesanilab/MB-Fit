import sys
from database import Database

def clean_database(settings_file, database_name):
    """
    Calls the database.clean() method on the given database. Sets all running calculations back to pending.
        
    """

    with Database(database_name) as database:

        database.clean()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python database_cleaner.py <settings.ini> <database_name>")
        exit(1)
    clean_database(sys.argv[1], sys.argv[2])
