"""
A module to figure out how to develop databases with sqlite3.
"""
import sqlite3
from datetime import datetime

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
        Molecules(mID int, mol blob, model text, package text, 
        energies blob, updated date)''')
    return cursor, connect

def query(cursor, table, obj_tuple, update=False):
    """
    Looks for an existing molecule configuration and do things with it
    Input: The cursor object, the name of the table, and the configuration,
           and a boolean for updating
    Output: The data that shall be returned if no inserting/updating
    """
    cursor.execute('''select * from {} where mol=?'''.format(table),
        (obj_tuple[1],))
    data = cursor.fetchone()

    if not data:
        insert(cursor, table, obj_tuple)
    elif update:
        update(cursor, table, obj_tuple)
    else:
        print("Retrieving current data")
        return data


def insert(cursor, table, obj_tuple):
    """
    Inserts an object tuple into the table
    Input: The cursor object, the name of the table, and the values
           required by the tuple itself
    """
    print("New configuration")
    cursor.execute('''insert into {} values{}'''.format(table, obj_tuple))

def update(cursor, table, obj_tuple):
    """
    Updates a certain entry in the database
    """
    print("User intends to update energy calculation")
    cursor.execute('''update {} set energies=?, updated=? where mol=?'''
        .format(table), (obj_tuple[2], datetime.now(), obj_tuple[1],))


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
insert(cursor,"Molecules",(1,0,0,0), True)
finalize(connect)
'''
