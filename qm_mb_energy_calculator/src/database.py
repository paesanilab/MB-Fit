"""
A module to figure out how to develop databases with sqlite3.
"""
import sqlite3
import re
import pickle
import sys

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

def query_id(cursor, table, id):
    """
    Looks for an existing molecule configuration and do things with it
    Input: The cursor object, the name of the table, and some arguments
    Output: The data that shall be returned, or None
    """

    cursor.execute('''select * from {} where ID=?'''.format(table), (id,))
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
        for char in key:
            if char not in invalid:
                newkey += char
            else:
                is_tuple = True
        if is_tuple:
            newkey = 'E' + newkey
        
        columns.append(newkey)
        entries.append(value)

    columns = tuple(columns)
    entries = tuple(entries)

    #print(process_string(entries))

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
    def __init__(self, file_name):
        self.connection = sqlite3.connect(file_name)
        self.cursor = self.connection.cursor()

        try:
            self.cursor.execute("PRAGMA table_info('schema_version')");
        except:
            raise ValueError("{} exists but is not a valid database file. \n Terminating database initialization.".format(database_name))
            sys.exit(1)

    '''
    Saves changes to the database
    '''
    def save(self):
        self.connection.commit()

    '''
    Closes the database
    '''
    def close(self):
        self.connection.close()

    '''
    Create the tables in the database
    '''
    def create(self):
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS Configs(ID text, config BLOB, natoms INT,
            nfrags INT)
            """)
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS Energies(ID TEXT, method TEXT, basis TEXT, cp TEXT, tag TEXT,
            E0 REAL, E1 REAL, E2 REAL, E01 REAL, E02 REAL, E12 REAL, E012 REAL, Enb BLOB)
            """)

    '''
    Add rows to the database
    '''
    def add_molecule(self, molecule, method, basis, cp, tag):
        hash_id = molecule.get_SHA1()

        self.cursor.execute("select ID from Configs where ID=?", (hash_id,))
        config_row = self.cursor.fetchone()
        if config_row is None:
            config = pickle.dumps(molecule).hex()
            self.cursor.execute("INSERT INTO Configs (ID, config, natoms, nfrags) VALUES ('{}', '{}', '{}', '{}')".format(hash_id, config, molecule.get_num_atoms(), molecule.get_num_fragments()))

        
        self.cursor.execute("select ID from Energies where ID=? AND method=? AND basis=? AND cp=? AND tag=?", (hash_id, method, basis, cp, tag))
        energies_row = self.cursor.fetchone()

        # create a new row in the Energies table if such a row did not already exist
        if energies_row is None:
            fragment_count = molecule.get_num_fragments();
            if fragment_count == 3:
                self.cursor.execute("INSERT INTO Energies (ID, method, basos, cp, tag, E0, E1, E2, E01, E02, E12, E012, Enb) VALUES ('{}', '{}', '{}', '{}', '{}', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None')".format(hash_id, method, basis, cp, tag))
            elif fragment_count == 2:
                self.cursor.execute("INSERT INTO Energies (ID, method, basis, cp, tag, E0, E1, E2, E01, E02, E12, E012, Enb) VALUES ('{}', '{}', '{}', '{}', '{}', 'None', 'None', 'N/A', 'None', 'N/A', 'N/A', 'N/A', 'None')".format(hash_id, method, basis, cp, tag))
            elif fragment_count == 1:
                self.cursor.execute("INSERT INTO Energies (ID, method, basis, cp, tag, E0, E1, E2, E01, E02, E12, E012, Enb) VALUES ('{}', '{}', '{}', '{}', '{}', 'None', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'None')".format(hash_id, method, basis, cp, tag))
            else:
                raise ValueError("Unsupported Number of fragments {}. Supported values are 1,2, and 3.".format(fragment_count))


    '''
    Gets one row with a mising energy
    '''
    def get_missing_energy(self):
        self.cursor.execute("SELECT * FROM Energies WHERE E0='None' OR E1='None' OR E2='None' OR E01='None' OR E12='None' OR E02='None' OR E012='None'")
        row = self.cursor.fetchone()
        if row is None:
            return None

        ID = row[0]
        method = row[1]
        basis = row[2]
        cp = row[3]
        tag = row[4]
        fragments = [0] if row[5] == "None" else [1] if row[6] == "None" else [2] if row[7] == "None" else [0,1] if row[8] == "None" else [1,2] if row[9] == "None" else [0,2] if row[10] == "None" else [0,1,2]
       
        self.cursor.execute("SELECT * FROM Configs WHERE ID=?", (ID,))
        molecule = pickle.loads(bytes.fromhex(self.cursor.fetchone()[1]))
        return Calculation(molecule, method, basis, cp, tag, fragments, None)
        

    '''
    Set an energy in the database
    '''
    def set_energy(self, calculation):
        entry_string = "E"
        for index in calculation.fragments:
            entry_string += str(index)
        self.cursor.execute("UPDATE Energies SET {}=? WHERE ID=? AND method=? AND basis=? AND cp=? AND tag=?".format(entry_string), (calculation.energy, calculation.molecule.get_SHA1(), calculation.method, calculation.basis, calculation.cp, calculation.tag))

    '''
    Gets list of all the Molecules with calculated energies
    TODO
    '''
    def get_complete_calculations(self):
        # get a list of all the rows with filled energies
        self.cursor.execute("SELECT * FROM Energies WHERE NOT (E0='None' OR E1='None' OR E2='None' OR E01='None' OR E12='None' OR E02='None' OR E012='None')")
        rows = self.cursor.fetchall()

        calculations = []

        for row in rows:
            ID = row[0]
            method = row[1]
            basis = row[2]
            cp = row[3]
            tag = row[4]
            fragments = [0] if row[5] == "None" else [1] if row[6] == "None" else [2] if row[7] == "None" else [0,1] if row[8] == "None" else [1,2] if row[9] == "None" else [0,2] if row[10] == "None" else [0,1,2]
       
            self.cursor.execute("SELECT * FROM Configs WHERE ID=?", (ID,))
            molecule = pickle.loads(bytes.fromhex(self.cursor.fetchone()[1]))

            calculations.append(Calculation(molecule, method, basis, cp, tag, fragments, energy))

        return calculations

class Calculation():
    def __init__(self, molecule, method, basis, cp, tag, fragments, energy):
        self.molecule = molecule
        self.method = method
        self.basis = basis
        self.cp = cp
        self.tag = tag
        self.fragments = fragments
        self.energy = energy
