"""
Contains the Database class, used to read and write to a database
"""
import sqlite3
import pickle
import sys

class Database():
    """
    Database class. Allows one to access a database and perform operations on it.
    """
    
    def __init__(self, file_name):
        """
        Initializer, sets up connection, and cursor
        """
        
        # connection is used to get the cursor, commit to the database, and close the database
        self.connection = sqlite3.connect(file_name)
        
        # the cursor is used to execute operations on the database
        self.cursor = self.connection.cursor()

        # this checks to make sure the file specified by file_name is a valid database file. The sqlite3 command will throw an error if the file is exists but is not a database
        try:
            self.cursor.execute("PRAGMA table_info('schema_version')");
        except:
            self.close()
            raise ValueError("{} exists but is not a valid database file. \n Terminating database initialization.".format(database_name)) from None

    def save(self):
        """
        Save any Changes to the database.
        Changes will not be perminent untill this function is called.
        """

        self.connection.commit()

    def close(self):
        """
        Close the database.
        Always close the database after you are done using it.
        """

        self.connection.close()

    def create(self):
        """
        Creates the required tables in the database, if they do not exists
        """
        
        # create the Configs table, which contains ID (SHA1 hash of the molecule), config (molecule data compressed as hexidecimal string), natoms (number of atoms), nfrags (number of fragments).
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS Configs(ID text, config BLOB, natoms INT,
            nfrags INT)
            """)

        # create the Energies table, which contains ID (SHA1 hash of the molecule), method, basis, cp (wheter counterpoise correction was used), tag (way to mark certain calculations), E0-E012 (nmer energies of the fragments in this molecule), and Enb (nbody energies)
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS Energies(ID TEXT, method TEXT, basis TEXT, cp TEXT, tag TEXT,
            E0 REAL, E1 REAL, E2 REAL, E01 REAL, E02 REAL, E12 REAL, E012 REAL, Enb BLOB)
            """)

    def add_molecule(self, molecule, method, basis, cp, tag):
        """
        Add a molecule to the database to have its energies calculated.
        """
        hash_id = molecule.get_SHA1()

        # find if a row exists for this molecule in the Configs table
        self.cursor.execute("select ID from Configs where ID=?", (hash_id,))
        config_row = self.cursor.fetchone()

        # if there is no row for this molecule in the Configs table, create one
        if config_row is None:
            config = pickle.dumps(molecule).hex() # Compress molecule data
            self.cursor.execute("INSERT INTO Configs (ID, config, natoms, nfrags) VALUES ('{}', '{}', '{}', '{}')".format(hash_id, config, molecule.get_num_atoms(), molecule.get_num_fragments()))

        # find if a row exists for this molecule with this specific set of method/basis/cp/tag in the Energies table
        self.cursor.execute("select ID from Energies where ID=? AND method=? AND basis=? AND cp=? AND tag=?", (hash_id, method, basis, cp, tag))
        energies_row = self.cursor.fetchone()

        # create a new row in the Energies table if such a row did not already exist
        if energies_row is None:
            fragment_count = molecule.get_num_fragments();

            # energies that do not apply to a molecule, due to having less than three fragments, are set to N/A, while all other energies are set to None
            if fragment_count == 3:
                self.cursor.execute("INSERT INTO Energies (ID, method, basis, cp, tag, E0, E1, E2, E01, E02, E12, E012, Enb) VALUES ('{}', '{}', '{}', '{}', '{}', 'None', 'None', 'None', 'None', 'None', 'None', 'None', 'None')".format(hash_id, method, basis, cp, tag))
            elif fragment_count == 2:
                self.cursor.execute("INSERT INTO Energies (ID, method, basis, cp, tag, E0, E1, E2, E01, E02, E12, E012, Enb) VALUES ('{}', '{}', '{}', '{}', '{}', 'None', 'None', 'N/A', 'None', 'N/A', 'N/A', 'N/A', 'None')".format(hash_id, method, basis, cp, tag))
            elif fragment_count == 1:
                self.cursor.execute("INSERT INTO Energies (ID, method, basis, cp, tag, E0, E1, E2, E01, E02, E12, E012, Enb) VALUES ('{}', '{}', '{}', '{}', '{}', 'None', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'None')".format(hash_id, method, basis, cp, tag))
            else:
                raise ValueError("Unsupported Number of fragments {}. Supported values are 1,2, and 3.".format(fragment_count))


    def get_missing_energy(self):
        """
        Returns a single Calculation object, which contains the info a user needs to calculate an energy missing in the table.

        The user should calculate the energy, set the calculation.energy field to the resultant energy, and then call set_energy(calculation)
        """

        # find a single row missing at least one energy
        self.cursor.execute("SELECT * FROM Energies WHERE E0='None' OR E1='None' OR E2='None' OR E01='None' OR E12='None' OR E02='None' OR E012='None'")
        row = self.cursor.fetchone()

        # if there is no such row, then the database is complete!
        if row is None:
            return None

        ID = row[0]
        method = row[1]
        basis = row[2]
        cp = row[3]
        tag = row[4]

        # this is a messy hardcoding, but essentially it gets an array of the fragments needed to calculate a missing energy
        fragments = [0] if row[5] == "None" else [1] if row[6] == "None" else [2] if row[7] == "None" else [0,1] if row[8] == "None" else [0,2] if row[9] == "None" else [1,2] if row[10] == "None" else [0,1,2]
       
        # retrieve and uncompress the molecule from the Configs table
        self.cursor.execute("SELECT * FROM Configs WHERE ID=?", (ID,))
        molecule = pickle.loads(bytes.fromhex(self.cursor.fetchone()[1]))

        return Calculation(molecule, method, basis, cp, tag, fragments, None)

    def missing_energies(self):
        """
        A generator to generate all the missing energies
        """
        while True:
            calculation = self.get_missing_energy()
            
            if calculation is None:
                break

            yield calculation
        

    def set_energy(self, calculation):
        """
        Sets the energy in the table of a certain calculation
        """

        # create a string representing the energy to set from the fragments used to calculate this energy
        entry_string = "E"
        for index in calculation.fragments:
            entry_string += str(index)

        # update the energy value in the corresponding row of the table
        self.cursor.execute("UPDATE Energies SET {}=? WHERE ID=? AND method=? AND basis=? AND cp=? AND tag=?".format(entry_string), (calculation.energy, calculation.molecule.get_SHA1(), calculation.method, calculation.basis, calculation.cp, calculation.tag))

# should probably be changed to a generator once we find a way to store the optimized geometry
    def get_complete_energies(self):
        """
        Returns a list of pairs of [molecule, energies] where energies is an array of the form [E0, E1, E2, E01, E12, E02, E012], where N/A energies are left out

        """

        # get a list of all the rows with completely filled energies
        self.cursor.execute("SELECT * FROM Energies WHERE NOT (E0='None' OR E1='None' OR E2='None' OR E01='None' OR E12='None' OR E02='None' OR E012='None')")
        rows = self.cursor.fetchall()

        molecule_energy_pairs = []

        for row in rows:
            # retrieve and uncompress the molecule from the corresponding Configs row
            self.cursor.execute("SELECT * FROM Configs WHERE ID=?", (row[0],))
            molecule = pickle.loads(bytes.fromhex(self.cursor.fetchone()[1]))
            num_frags = molecule.get_num_fragments()
            
            energies = []

            # fill the energies array with the correct values based on the number of fragments in this molecule. Monomers will have just 1 energy, dimers 3, and trimers 7
            if num_frags == 3:
                for i in range(5, 12):
                    energies.append(row[i])
            elif num_frags == 2:
                for i in [5, 6, 8]:
                    energies.append(row[i])
            elif num_frags == 1:
                energies.append(row[5])
            else:
                raise ValueError("Unsupported Number of fragments {}. Supported values are 1,2, and 3.".format(fragment_count))

            # add this [molecule, energies] pair to the list of pairs
            molecule_energy_pairs.append([molecule, energies])

        return molecule_energy_pairs

    def get_comparison_energies(self, method1, method2, basis1, basis2, cp1, cp2, energy):
        """
        Allows the user to search the database for two different sets of data corresponding to each set of method/basis/cp and
        Return an array of pairs of energies of the same molecule calculated in each way
        """
         
        # get rows pertaining to category 1
        self.cursor.execute("SELECT ID FROM Energies WHERE method=? AND basis=? AND cp=?", (method1, basis1, cp1))
        category1_IDs = self.cursor.fetchall()
        self.cursor.execute("SELECT {} FROM Energies WHERE method=? AND basis=? AND cp=?".format(energy), (method1, basis1, cp1))
        category1_energies = self.cursor.fetchall()

        # get rows pertaining to category 2
        self.cursor.execute("SELECT ID FROM Energies WHERE method=? AND basis=? AND cp=?", (method2, basis2, cp2))
        category2_IDs = self.cursor.fetchall()
        self.cursor.execute("SELECT {} FROM Energies WHERE method=? AND basis=? AND cp=?".format(energy), (method2, basis2, cp2))
        category2_energies = self.cursor.fetchall()

        # make list of all ids
        IDs = category1_IDs.copy()
        for ID in category2_IDs:
            if ID not in category1_IDs:
                IDs.append(ID)
        
        energy_pairs = []
       
        # fill energy pairs with pairs of energies corresponding to the same molecule in each category
        for ID in IDs:
            try:
                category1_energy = category1_energies[category1_IDs.index(ID)][0]
                category2_energy = category2_energies[category2_IDs.index(ID)][0]
            
                energy_pairs.append([category1_energy, category2_energy])
            except ValueError:
                pass

        return energy_pairs;


class Calculation():
    """
    Contains all the information the user needs to be able to make a calculation
    """
    def __init__(self, molecule, method, basis, cp, tag, fragments, energy):
        self.molecule = molecule
        self.method = method
        self.basis = basis
        self.cp = cp
        self.tag = tag
        self.fragments = fragments
        self.energy = energy
