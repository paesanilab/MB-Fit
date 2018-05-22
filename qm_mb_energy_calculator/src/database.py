"""
A module to figure out how to develop databases with sqlite3.
"""
import sqlite3
from datetime import datetime
import re
from exceptions import MoleculeIDUnsetException
from exceptions import RowOperationBeforeInitException

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

    columns = []
    entries = []
    invalid = '(), '
    is_tuple = False
    
    for key, value in kwargs.items():

        # Build new string to conform with sqlite3
        newkey = ''
        print("Before: "+key)
        for char in key:
            if char not in invalid:
                newkey += char
            else:
                is_tuple = True
        if is_tuple:
            newkey = 'E' + newkey
        
        # Convert value to conform with sqlite3
        if isinstance(value, (list,)):
            value = tuple(value)

        columns.append(newkey)
        entries.append(value)

        print("After: "+newkey)

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

    cursor.execute('''update {} set {} where ID={}'''.format(table, update_input, hashID))

def finalize(connection):
    """
    Ends the current session by saving and closing the database.
    Input: The connection to the database
    """
    connection.commit()
    connection.close()

"""
Data base class
"""
class Database():
    '''
    Initialize the database
    '''
    def __init__(self, name):
        self.name = name
        self.cursor, self.connection = __init__(name)
        self.molecule_id = None

    '''
    Close the database
    '''
    def close(self):
        finalize(self.connection)

    '''
    Update the molecule ID, which dictates what molecule we are working with.
    Molecule ID MUST be set before any oeprations can be performed.
    '''
    def set_molecule_id(self, molecule_id):
        self.molecule_id = molecule_id

    '''
    Clears the molecule ID,
    '''
    def clear_molecule_id(self):
        self.molecule_id = None

    '''
    Checks if molecule ID is set
    '''
    def has_id(self):
        return id is not None

    '''
    Initializes the molecule in the database by creating tables and inserting
    basic information into the table

    MUST be called before a molecule's data can be added from the database
    '''
    def init_molecule(self, config, natoms, nfrags, tag_in, model_in, cp_in):
        # check if molecule id is initialized
        if not self.has_id():
            raise MoleculeIDUnsetException("init_molecule()", self.name)

        # check if molecule already has rows in the table
        self.cursor.execute("select ID from Configs where ID=?", (self.molecule_id,))
        config_row = self.cursor.fetchone()
        self.cursor.execute("select ID from Energies where ID=?", (self.molecule_id,))
        energies_row = self.cursor.fetchone()
        
        # create a new row in the Configs table if such a row did not already exist
        if config_row is None: 
            insert(self.cursor, "Configs", ID=self.molecule_id, config="temp", # TODO: buffer doesn't work, "temp" should be buffer(compressed_mol)
                natom=natoms, nfrags=nfrags, tag=tag_in)

        # create a new row in the Energies table if such a row did not already exist
        if energies_row is None:
            insert(self.cursor, "Energies", ID=self.molecule_id, model=model_in, cp=cp_in, E0="UNSET", E1="UNSET",
                E2="UNSET", E01="UNSET", E02="UNSET", E12="UNSET", E012="UNSET", Enb="UNSET")  

    '''
    Checks if a certain molecule's n-mer energy is already in the database
    '''
    def has_nmer_energy(self, nmer):
        # check if molecule id is initialized
        if not self.has_id():
            raise MoleculeIDUnsetException("has_nmer_energy", self.name)

        # makes nmer string out of combination passed in, at the end, format should be "E1", "E23", etc.
        nmer_string = "E"
        for index in nmer:
            nmer_string += str(index)

        # fetch nmer energy from database
        self.cursor.execute("select {} from Energies where ID=?".format(nmer_string), (self.molecule_id,))
        energy = self.cursor.fetchone()

        # check if row was never initialized
        if energy is None:
            raise RowOperationBeforeInitException("has_nmer_energy", self.name)
        
        # return whether energy of nmer is set
        return energy[0] != "UNSET"
        

    '''
    Set's a certain molecule's n-mer energy, or resets it if already set
    '''
    def set_nmer_energy(self, nmer, energy):
        # check if molecle id is initialized
        if not self.has_id():
            raise MoleculeIDUnsetException("set_nmer_energy", self.name)

        # makes nmer string out of combination passed in, at the end, format should be "E1", "E23", etc.
        nmer_string = "E"
        for index in nmer:
            nmer_string += str(index)

        # first check if energy row exists in database
        self.cursor.execute("select {} from Energies where ID=?".format(nmer_string), (self.molecule_id,))
        energy_in_table = self.cursor.fetchone()
        if energy_in_table is None:
            raise RowbaseOperationBeforeInitException("set_nmer_energy", self.name)

        # if row exists, then go ahead and update table entry
        self.cursor.execute("update Energies set {}=? where ID=?".format(nmer_string), (str(energy), self.molecule_id))

    '''
    Get's a certain molecule's n-mer energy
    '''
    def get_nmer_energy(self, nmer):
        # check if molecle id is initialized
        if not self.has_id():
            raise MoleculeIDUnsetException("get_nmer_energy", self.name)

        # makes nmer string out of combination passed in, at the end, format should be "E1", "E23", etc.
        nmer_string = "E"
        for index in nmer:
            nmer_string += str(index)

        # first check if energy row exists in database
        self.cursor.execute("select {} from Energies where ID=?".format(nmer_string), (self.molecule_id,))
        energy_in_table = self.cursor.fetchone()
        if energy_in_table is None:
            raise RowbaseOperationBeforeInitException("set_nmer_energy", self.name)

        # TODO: some sort of check to see if energy is unset, if energy is unset, throw exception

        # if row exists, return energy
        return energy_in_table[0]
        

# All lines below functional
'''
cursor, connect = __init__("a.db")
insert(cursor,"Molecules",(1,0,0,0,0))
print(type(query(cursor,"Molecules",(1,0,0,0,0)))) # fetchone returns a tuple
finalize(connect)
'''
