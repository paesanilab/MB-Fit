"""
A module to figure out how to develop databases with sqlite3.
"""
import sqlite3

connect = sqlite3.connect("a.db")
c = connect.cursor()
#c.execute('''drop table if exists Molecules''')
c.execute('''create table if not exists 
    Molecules(mID int, mol blob, energies blob, updated date)''')
c.execute("insert into Molecules values(2,0,0,0)")
connect.commit()
connect.close()

def __init__(database_name):
    """
    Initializing the database for the user to work on.
    Input: Name of the database desired.
    Output: A cursor object to perform query commands, and the connection to it.
    """
    # Attach .db to the end if there does not exist any
    connect = sqlite3.connect(database_name)
    cursor = connect.cursor()
    cursor.execute('''create table if not exists 
    Molecules(mID int, mol blob, energies blob, updated date)''')
    return cursor, connect

def insert(cursor, table, obj_tuple):
    """
    Inserts an object tuple into the table
    Input: The cursor object, the name of the table, and the values
           required by the tuple itself
    """
    cursor.execute('''insert into {} values{}'''.format(table, obj_tuple))

def finalize(connection)
    """
    Ends the current session by saving and closing the database.
    Input: The connection to the database
    """
    connection.commit()
    connection.close()
