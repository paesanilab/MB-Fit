import sys

from mbfit.database import Database

if(len(sys.argv) != 2):
    print("Usage:")
    print("{} <database config>".format(sys.argv[0]))
    exit(1)

database_config = sys.argv[1]

with Database(database_config) as database:
    database.create()
