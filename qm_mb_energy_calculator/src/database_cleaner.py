import sys
from database import Database

def clean_database(settings_file, database_name):
    """
    Calls the database.clean() method on the given database. Sets all running calculations back to pending.
        
    """
    # add .db to database name if it doesn't already end in .db 
    if database_name[-3:] != ".db":
        print("Database name \"{}\" does not end in database suffix \".db\". Automatically adding \".db\" to end of database name.".format(database_name))
        database_name += ".db"

    with Database(database_name) as database:

        database.clean()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python database_cleaner.py <settings.ini> <database_name>")
        exit(1)
    clean_database(sys.argv[1], sys.argv[2])
