"""
Contains the Database class, used to read and write to a database
"""
import sqlite3
import pickle
import sys
from molecule import Atom, Fragment, Molecule
import itertools
import datetime

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
        
        # create the Models table 
        self.cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS Models(method TEXT, basis TEXT, cp INT)
            """
        )

        # create the Molecules table
        self.cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS Molecules(hash TEXT)
            """
        )

        # create the Fragments table
        self.cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS Fragments(molecule_id INT, charge INT, spin INT)
            """
        )

        # create the Atoms table
        self.cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS Atoms(fragment_id INT, symbol TEXT, x REAL, y REAL, z REAL)
            """
        )

        # create the Calculations table
        self.cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS Calculations(molecule_id INT, model_id INT, tag TEXT, optimized INT)
            """
        )

        # create the Energies table
        self.cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS Energies(calculation_id INT, job_id INT, energy_index INT, energy REAL)
            """
        )

        # create the Jobs table
        self.cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS Jobs(status TEXT, log_file TEXT, start_date TEXT, end_date TEXT)
            """
        )

    def add_calculation(self, molecule, method, basis, cp, tag, optimized):
        """
        Add a molecule to the database to have its energies calculated.
        """

        # check if this model is not already in Models table
        if not self.cursor.execute("SELECT EXISTS(SELECT * FROM Models WHERE method=? AND basis=? AND cp=?)", (method, basis, cp)).fetchone()[0]:

            # create entry in Models table
            self.cursor.execute("INSERT INTO Models (method, basis, cp) VALUES (?, ?, ?)", (method, basis, cp))
            
            # get id of this model
            model_id = self.cursor.lastrowid

        else:

            # get id of this model
            model_id = self.cursor.execute("SELECT ROWID FROM Models WHERE method=? AND basis=? AND cp=?", (method, basis, cp)).fetchone()[0]

        # get the SHA1 hash of this molecule
        molecule_hash = molecule.get_SHA1()

        # check if this molecule is not already in Molecules table
        if not self.cursor.execute("SELECT EXISTS(SELECT * FROM Molecules WHERE hash=?)", (molecule_hash,)).fetchone()[0]:

            # create entry in Molecules table
            self.cursor.execute("INSERT INTO Molecules (hash) VALUES (?)", (molecule_hash,))
        
            # get id of this molecule
            molecule_id = self.cursor.lastrowid

            # insert molecule's fragments into the table
            for fragment in molecule.get_fragments():
                self.cursor.execute("INSERT INTO Fragments (molecule_id, charge, spin) VALUES (?, ?, ?)", (molecule_id, fragment.get_charge(), fragment.get_spin_multiplicity()))
                
                # get id of this fragment
                fragment_id = self.cursor.lastrowid

                # insert fragment's atoms into the table
                for atom in fragment.get_atoms():
                   self.cursor.execute("INSERT INTO Atoms (fragment_id, symbol, x, y, z) VALUES (?, ?, ?, ?, ?)", (fragment_id, atom.get_name(), atom.get_x(), atom.get_y(), atom.get_z())) 
                    
        else:
            
            # get id of this molecule
            molecule_id = self.cursor.execute("SELECT ROWID FROM Molecules WHERE hash=?", (molecule_hash,)).fetchone()[0]

        # check if the calculation is not already in the Calculations table

        if not self.cursor.execute("SELECT EXISTS(SELECT * FROM Calculations WHERE molecule_id=? AND model_id=? AND tag=? AND optimized=?)", (molecule_id, model_id, tag, optimized)).fetchone()[0]:
            
            # create entry in Calculations table
            self.cursor.execute("INSERT INTO Calculations (molecule_id, model_id, tag, optimized) VALUES (?, ?, ?, ?)", (molecule_id, model_id, tag, optimized))

            # get id of this calculation
            calculation_id = self.cursor.lastrowid

            # add rows to the Energies table
            for energy_index in range(number_of_energies(molecule.get_num_fragments())):
                # create a job for this energy
                self.cursor.execute("INSERT INTO Jobs (status) VALUES (?)", ("pending",))

                # get the id of this job
                job_id = self.cursor.lastrowid
        
                # insert row into Energies table for this energy
                self.cursor.execute("INSERT INTO Energies (calculation_id, energy_index, job_id) VALUES (?, ?, ?)", (calculation_id, energy_index, job_id))

    def get_missing_energy(self):
        """
        Returns a single Calculation object, which contains the info a user needs to calculate an energy missing in the table.

        The user should calculate the energy, set the calculation.energy field to the resultant energy, and then call set_energy(calculation)
        """

        # retrieve a pending job from the Jobs table
        try:
            job_id = self.cursor.execute("SELECT ROWID FROM Jobs WHERE status=?", ("pending",)).fetchone()[0]
        except TypeError:
            # this error will be thrown if there are no more energies to calculate, because None[0] throws a TypeError
            return None
        
        # retrieve the calculation and energy to be calculated for this job from the Energies table
        calculation_id, energy_index = self.cursor.execute("SELECT calculation_id, energy_index FROM Energies WHERE job_id=?", (job_id,)).fetchone()

        # retrieve the molecule and model from the Calculations table
        molecule_id, model_id, tag, optimized = self.cursor.execute("SELECT molecule_id, model_id, tag, optimized FROM Calculations WHERE ROWID=?", (calculation_id,)).fetchone()

        # retrieve the method, basis, and cp for this model
        method, basis, cp = self.cursor.execute("SELECT method, basis, cp FROM Models WHERE ROWID=?", (model_id,)).fetchone()

        # Reconstruct the molecule from the information in this database
        molecule = Molecule()
        
        # loop over all rows in the Fragments table that correspond to this molecule
        for fragment_id, charge, spin in self.cursor.execute("SELECT ROWID, charge, spin FROM Fragments WHERE molecule_id=?", (molecule_id,)).fetchall():
            fragment = Fragment(charge, spin)

            # loop over all rows in the Atoms table that correspond to this fragment
            for symbol, x, y, z in self.cursor.execute("SELECT symbol, x, y, z FROM Atoms WHERE fragment_id=?", (fragment_id,)).fetchall():
                fragment.add_atom(Atom(symbol, x, y, z))

            molecule.add_fragment(fragment)

        # get the indicies of the fragments to include in this calculation
        fragment_indicies = energy_index_to_fragment_indicies(energy_index, molecule.get_num_fragments())

        # update the job with its start date and running status
        self.cursor.execute("UPDATE Jobs SET status=?, start_date=? WHERE ROWID=?", ("running", datetime.datetime.today().strftime('%Y/%m/%d'), job_id))
        
        return Calculation(molecule, method, basis, True if cp == 1 else False, tag, True if optimized == 1 else False, fragment_indicies)

    def missing_energies(self):
        """
        A generator to generate all the missing energies
        """
        while True:
            calculation = self.get_missing_energy()
            
            if calculation is None:
                raise StopIteration

            yield calculation
        

    def set_energy(self, calculation, energy, log_file):
        """
        Sets the energy in the table of a certain calculation
        """

        # find the id of the model from the Models table
        try:
            model_id = self.cursor.execute("SELECT ROWID FROM Models WHERE method=? AND basis=? AND cp=?", (calculation.method, calculation.basis, calculation.cp)).fetchone()[0]
        except TypeError:
            # this error gets raised if this model is not in the Models table
            raise ValueError("Tried to set an energy that does not exist in the database") from None

        # get the SHA1 hash of this molecule
        molecule_hash = calculation.molecule.get_SHA1()

        # find the id of the molecule from the Molecules table
        try:
            molecule_id = self.cursor.execute("SELECT ROWID FROM Molecules WHERE hash=?", (molecule_hash,)).fetchone()[0]
        except TypeError:
            # this error gets raised if this molecule is not in the Molecules table
            raise ValueError("Tried to set an energy that does not exist in the database") from None        

        # find the id of the calculation from the Calculations table
        try:
            calculation_id = self.cursor.execute("SELECT ROWID FROM Calculations WHERE molecule_id=? AND model_id=? AND tag=? AND optimized=?", (molecule_id, model_id, calculation.tag, calculation.optimized)).fetchone()[0]
        except TyperError:
            # this error gets raised if this calculation is not in the Calculations table
            raise ValueError("Tried to set an energy that does not exist in the database")

        energy_index = fragment_indicies_to_energy_index(calculation.fragments, calculation.molecule.get_num_fragments())

        # update the row in the Energies table corresponding to this energy
        self.cursor.execute("UPDATE Energies SET energy=? WHERE calculation_id=? AND energy_index=?", (energy, calculation_id, energy_index))

        # get the job corresponding to this energy
        job_id = self.cursor.execute("SELECT job_id FROM Energies WHERE calculation_id=? AND energy_index=?", (calculation_id, energy_index)).fetchone()[0]

        # update the information about this job
        self.cursor.execute("UPDATE Jobs SET status=?, log_file=?, end_date=? WHERE ROWID=?", ("completed", log_file, datetime.datetime.today().strftime('%Y/%m/%d'), job_id))

    def get_complete_energies(self):
        """
        Returns a list of pairs of [molecule, energies] where energies is an array of the form [E0, ...]

        Includes any optimized geometries        
        """

        yield from self.get_energies("%", "%", "%", "%")

    def get_energies(self, method, basis, cp, tag):
        """
        Returns a list of pairs of [molecule, energies] where energies is an array of the form [E0, ...] of molecules in the database with the given
        method, basis, cp and tag.

        Includes any optimized geometries.

        % can be used as a wildcard to stand in for any method, basis, cp, or tag. 
        """

        # get a list of all calculations that have the appropriate method, basis, and cp
        calculation_ids = [elements[0] for elements in self.cursor.execute("SELECT ROWID FROM Calculations WHERE model_id=(SELECT ROWID FROM Models WHERE method LIKE ? AND basis LIKE ? AND cp LIKE ? AND tag LIKE ?)", (method, basis, cp, tag)).fetchall()]

        for calculation_id in calculation_ids:
            
            # get the molecule id corresponding to this calculation
            molecule_id = self.cursor.execute("SELECT molecule_id FROM Calculations WHERE ROWID=?", (calculation_id,)).fetchone()[0]

            # Reconstruct the molecule from the information in this database
            molecule = Molecule()
            
            # loop over all rows in the Fragments table that correspond to this molecule
            for fragment_id, charge, spin in self.cursor.execute("SELECT ROWID, charge, spin FROM Fragments WHERE molecule_id=?", (molecule_id,)).fetchall():
                fragment = Fragment(charge, spin)

                # loop over all rows in the Atoms table that correspond to this fragment
                for symbol, x, y, z in self.cursor.execute("SELECT symbol, x, y, z FROM Atoms WHERE fragment_id=?", (fragment_id,)).fetchall():
                    fragment.add_atom(Atom(symbol, x, y, z))

                molecule.add_fragment(fragment)
            
            # get the energies corresponding to this calculation
            energies = [energy_index_value_pair[1] for energy_index_value_pair in sorted(self.cursor.execute("SELECT energy_index, energy FROM Energies WHERE calculation_id=?", (calculation_id,)).fetchall())]
    
            yield molecule, energies

    def get_complete_optimized_energies(self):
        """
        Returns a list of pairs of [molecule, energies] where energies is an array of the form [E0, ...] for optimized geometries.
        """
        yield from self.get_optimized_energies("%", "%", "%", "%")

    def get_optimized_energies(self, method, basis, cp, tag):
        """
        Returns a list of pairs of [molecule, energies] where energies is an array of the form [E0, ...] of molecules in the database with the given
        method, basis, cp and tag and optimized geometries

        % can be used as a wildcard to stand in for any method, basis, cp, or tag. 
        """

        # get a list of all calculations that have the appropriate method, basis, and cp
        calculation_ids = [elements[0] for elements in self.cursor.execute("SELECT ROWID FROM Calculations WHERE model_id=(SELECT ROWID FROM Models WHERE method LIKE ? AND basis LIKE ? AND cp LIKE ? AND tag LIKE ? AND optimized=?)", (method, basis, cp, tag, 1)).fetchall()]

        for calculation_id in calculation_ids:
            
            # get the molecule id corresponding to this calculation
            molecule_id = self.cursor.execute("SELECT molecule_id FROM Calculations WHERE ROWID=?", (calculation_id,)).fetchone()[0]

            # Reconstruct the molecule from the information in this database
            molecule = Molecule()
            
            # loop over all rows in the Fragments table that correspond to this molecule
            for fragment_id, charge, spin in self.cursor.execute("SELECT ROWID, charge, spin FROM Fragments WHERE molecule_id=?", (molecule_id,)).fetchall():
                fragment = Fragment(charge, spin)

                # loop over all rows in the Atoms table that correspond to this fragment
                for symbol, x, y, z in self.cursor.execute("SELECT symbol, x, y, z FROM Atoms WHERE fragment_id=?", (fragment_id,)).fetchall():
                    fragment.add_atom(Atom(symbol, x, y, z))

                molecule.add_fragment(fragment)
            
            # get the energies corresponding to this calculation
            energies = [energy_index_value_pair[1] for energy_index_value_pair in sorted(self.cursor.execute("SELECT energy_index, energy FROM Energies WHERE calculation_id=?", (calculation_id,)).fetchall())]
    
            yield molecule, energies

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
    def __init__(self, molecule, method, basis, cp, tag, optimized, fragments):
        self.molecule = molecule
        self.method = method
        self.basis = basis
        self.cp = cp
        self.tag = tag
        self.optimized = optimized
        self.fragments = fragments

def number_of_energies(number_of_fragments):
    """
    Returns the number of nmer energies that a molecule with number_of_fragments fragments will have
    1 -> 1
    2 -> 3
    3 -> 7
    4 -> 15

    Equal to the size of the power set of a set of size number_of_fragments minus the empty set
    """

    return 2 ** number_of_fragments - 1

# THESE TWO METHODS COULD PROBABLY BE BETTER, BUT THEY WORK

def energy_index_to_fragment_indicies(energy_index, number_of_fragments):
    """
    Returns an array of fragment indicies that should be included in a calculation with energy_index index in a molecule with number_of_fragments fragments
    """

    combinations = [inner for outer in (itertools.combinations(range(number_of_fragments), combination_size) for combination_size in range(1, number_of_fragments + 1)) for inner in outer]

    return combinations[energy_index]

def fragment_indicies_to_energy_index(fragment_indicies, number_of_fragments):
    """
    Returns the energy index corresponding to an energy using fragments fragment_indicies and number_of_fragments fragments
    """

    combinations = [inner for outer in (itertools.combinations(range(number_of_fragments), combination_size) for combination_size in range(1, number_of_fragments + 1)) for inner in outer]

    return combinations.index(fragment_indicies)
