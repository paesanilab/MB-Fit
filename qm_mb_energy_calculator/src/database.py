"""
A module to figure out how to develop databases with sqlite3.
"""
import sqlite3
from datetime import datetime

# Much of the sqlite3 commands are hard-coded, but we do not need to make
# things generic here (yet)

def __init__(database_name):
    """
    Initializing the database for the user to work on.
    Input: Name of the database desired.
    Output: A cursor object to perform query commands, and the connection to it.
    """
    # Attach .db to the end if there does not exist any
    if database_name.find(".db") == -1:
        database_name += ".db"

    connect = sqlite3.connect(database_name)
    cursor = connect.cursor()
    cursor.execute('''create table if not exists 
        Configs(ID text, config blob, natom int, nfrags int, tag text)''')
    cursor.execute('''create table if not exists 
        Energies(ID text, model text, cp int, E0 real, E1 real, E2
        real, E01 real, E02 real, E12 real, E012 real, 
        Enb blob)''')

    return cursor, connect

def query(cursor, table, **kwargs):
    """
    Looks for an existing molecule configuration and do things with it
    Input: The cursor object, the name of the table, and some arguments
    Output: The data that shall be returned, or None
    """
    query_input=""
    
    for key, value in kwargs.items():
        query_input += "{}={} and ".format(key, value)
    query_input = query_input[:-4]
      
    cursor.execute('''select * from {} where {}'''.format(table, query_input))
    return cursor.fetchone()



def insert(cursor, table, **kwargs):
    """
    Inserts an object tuple into the table
    Input: The cursor object, the name of the table, and the values
           required by the tuple itself
    """

    columns=[]
    entries=[]
    
    for key, value in kwargs.items():

        print(key)

        # Special string processing to remove special characters
        key = key.translate(key.maketrans("", "", "(),"))
        
        columns.append(key)
        entries.append(value)

    columns = tuple(columns)
    entries = tuple(entries)

    print(columns)
    print(entries)

    cursor.execute('''insert into {} {} values {}'''.format(table, 
        columns,entries))

def update(cursor, table, hashID, **kwargs):
    """
    Updates a certain entry in the database
    Input: A pointer to the database, the table name, the object hash,
           and the arguments to update
    """
    update_input = ""
    
    for key, value in kwargs.items():
        update_input += "{}={},".format(key, value)

    update_input = update_input[:-1]

    cursor.execute('''update {} set {} where ID={}'''
        .format(table, update_input, hashID))

def finalize(connection):
    """
    Ends the current session by saving and closing the database.
    Input: The connection to the database
    """
    connection.commit()
    connection.close()

# All lines below functional
'''
cursor, connect = __init__("a.db")
insert(cursor,"Molecules",(1,0,0,0,0))
print(type(query(cursor,"Molecules",(1,0,0,0,0)))) # fetchone returns a tuple
finalize(connect)
'''
