"""
Contains the Database class, used to read and write to a database
"""
import itertools, datetime
import sqlite3
from potential_fitting.molecule import Atom, Fragment, Molecule
from potential_fitting.exceptions import InconsistentDatabaseError, InvalidValueError

class Database():
    """
    Database class. Allows one to access a database and perform operations on it.
    """
    
    def __init__(self, file_name):
        """
        Initializer, sets up connection, and cursor

        Args:
            file_name   - path to the database file

        Returns:
            a new Database object
        """

        if file_name[-3:] != ".db":
            print("Automatically ending '.db' suffix to database name {}.".format(file_name))
            file_name += ".db"

        self.file_name = file_name
        
        # connection is used to get the cursor, commit to the database, and close the database
        self.connection = sqlite3.connect(file_name)
        
        # the cursor is used to execute operations on the database
        self.cursor = self.connection.cursor()

        # this checks to make sure the file specified by file_name is a valid database file. The sqlite3 command will throw an error if the file is exists but is not a database
        try:
            self.cursor.execute("PRAGMA table_info('schema_version')");
        except:
            self.close()
            raise InconsistentDatabaseError(file_name, "File exists but is not a valid database file.") from None
    

    # the __enter__() and __exit__() methods define a database as a context manager, meaning you can use with ... as ... syntax on it
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.save()
        self.close()
        
        # returning false lets the context manager know that no exceptions were handled in the __exit__() method
        return False

    def save(self):
        """
        Save any Changes to the database.
        Changes will not be perminent untill this function is called.

        Args:
            None
        
        Returns:
            None
        """

        self.connection.commit()

    def close(self):
        """
        Close the database.
        Always close the database after you are done using it.

        Args:
            None

        Returns:
            None
        """

        self.connection.close()

    def create(self):
        """
        Creates the all tables in the database, if they do not exists

        Args:
            None

        Returns:
            None
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
            CREATE TABLE IF NOT EXISTS Molecules(name TEXT, hash TEXT)
            """
        )

        # create the Fragments table
        self.cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS Fragments(molecule_id INT, name TEXT, charge INT, spin INT)
            """
        )

        # create the Atoms table
        self.cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS Atoms(fragment_id INT, symbol TEXT, symmetry_class TEST, x REAL, y REAL, z REAL)
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
        Add a calculation to the database.

        Args:
            molecule - the molecule to calculate the nmer energies of
            method  - the method to use to calculate the energies
            basis   - the basis to use to calculate the energies
            cp      - whether to use counterpoise correction for this calculation
            tag     - a special tag to be able to identify this calculation at a later date
            optimized - whether this calculation uses an optimized geometry

        Returns:
            None
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
        if not self.cursor.execute("SELECT EXISTS(SELECT * FROM Molecules WHERE name=? AND hash=?)", (molecule.get_name(), molecule_hash)).fetchone()[0]:

            # create entry in Molecules table
            self.cursor.execute("INSERT INTO Molecules (name, hash) VALUES (?, ?)", (molecule.get_name(), molecule_hash))
        
            # get id of this molecule
            molecule_id = self.cursor.lastrowid

            # insert molecule's fragments into the table
            for fragment in molecule.get_fragments():
                self.cursor.execute("INSERT INTO Fragments (molecule_id, name, charge, spin) VALUES (?, ?, ?, ?)", (molecule_id, fragment.get_name(), fragment.get_charge(), fragment.get_spin_multiplicity()))
                
                # get id of this fragment
                fragment_id = self.cursor.lastrowid

                # insert fragment's atoms into the table
                for atom in fragment.get_atoms():
                   self.cursor.execute("INSERT INTO Atoms (fragment_id, symbol, symmetry_class, x, y, z) VALUES (?, ?, ?, ?, ?, ?)", (fragment_id, atom.get_name(), atom.get_symmetry_class(), atom.get_x(), atom.get_y(), atom.get_z())) 
                    
        else:
            
            # get id of this molecule
            molecule_id = self.cursor.execute("SELECT ROWID FROM Molecules WHERE name=? AND hash=?", (molecule.get_name(), molecule_hash)).fetchone()[0]

        # check if the calculation is not already in the Calculations table

        if not self.cursor.execute("SELECT EXISTS(SELECT * FROM Calculations WHERE molecule_id=? AND model_id=? AND tag=? AND optimized=?)", (molecule_id, model_id, tag, optimized)).fetchone()[0]:
            
            # create entry in Calculations table
            self.cursor.execute("INSERT INTO Calculations (molecule_id, model_id, tag, optimized) VALUES (?, ?, ?, ?)", (molecule_id, model_id, tag, optimized))

            # get id of this calculation
            calculation_id = self.cursor.lastrowid

            # add rows to the Energies table
            for energy_index in range(number_of_energies(molecule.get_num_fragments(), cp)):
                # create a job for this energy
                self.cursor.execute("INSERT INTO Jobs (status) VALUES (?)", ("pending",))

                # get the id of this job
                job_id = self.cursor.lastrowid
        
                # insert row into Energies table for this energy
                self.cursor.execute("INSERT INTO Energies (calculation_id, energy_index, job_id) VALUES (?, ?, ?)", (calculation_id, energy_index, job_id))

    def get_missing_energy(self):
        """
        Returns a single Job object, which contains the info a user needs to calculate an energy missing in the table.

        The user should calculate the energy, then call either set_energy() or set_failed()

        Args:
            None

        Returns:
            A Job object describing the calculation to be performed
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
        molecule = self.get_molecule(molecule_id)

        # get the indicies of the fragments to include in this calculation
        fragment_indicies = energy_index_to_fragment_indicies(energy_index, molecule.get_num_fragments(), True if cp == 1 else False)

        # update the job with its start date and running status
        self.cursor.execute("UPDATE Jobs SET status=?, start_date=? WHERE ROWID=?", ("running", datetime.datetime.today().strftime('%Y/%m/%d'), job_id))
        
        return Job(molecule, method, basis, True if cp == 1 and energy_index < number_of_energies(molecule.get_num_fragments(), False) else False, fragment_indicies, job_id)

    def missing_energies(self):
        """
        A generator to generate Jobs for all the missing energies

        Args:
            None

        Returns:
            A generator which generates Job objects until one has been generated for every pending job in the database
        """
        while True:
            calculation = self.get_missing_energy()
            
            if calculation is None:
                return

            yield calculation
        

    def set_energy(self, job_id, energy, log_file):
        """
        Sets the energy in the table of a certain calculation

        Args:
            job_id  - the id of the job associated with this energy, included in the Job object retrieved by get_missing_energy()
            energy  - the calculated energy
            log_file - path to the log file for the job

        Returns:
            None
        """
        
        if self.cursor.execute("SELECT status FROM Jobs WHERE ROWID=?", (job_id,)).fetchone()[0] != "running":
            print("This job is not set to running, but tried to have its energy set. Most likely, this job was dispatched a long time ago, and the database has since given up on ever getting a response. Energy will be placed into the database anyways, this could create an issue as the currently running job may later update this energy.") # Unsure how to handle this? Throw Error? Allow energy into database or not?


        # update the row in the Energies table corresponding to this energy
        self.cursor.execute("UPDATE Energies SET energy=? WHERE job_id=?", (energy, job_id))

        # update the information about this job
        self.cursor.execute("UPDATE Jobs SET status=?, log_file=?, end_date=? WHERE ROWID=?", ("completed", log_file, datetime.datetime.today().strftime('%Y/%m/%d'), job_id))

    def set_failed(self, job_id, status, log_path):
        """
        Sets the status of a job to a specific status

        Args:
            job_id  - the id of the job to set the status of, included in the Job object retrieved by get_missing_energy()
            status  - the new status of this job. Valid statuses are pending, running, completed, failed
            log_file - path to the log file for the job

        Returns:
            None
        """
        if status not in ["pending", "running", "completed", "failed"]:
            raise InvalidValueError("status", status, "one of 'pending', 'running', 'completed', or 'failed'")

        self.cursor.execute("UPDATE Jobs SET status=?, log_file=? WHERE ROWID=?", (status, log_path, job_id))

    def get_complete_energies(self, molecule_name, optimized = False):
        """
        Returns a generator of pairs of [molecule, energies] where energies is an array of the form [E0, ...]

        Includes any optimized geometries

        Args:
            None

        Returns:
            a generator of [molecule, [E0, E1, E2, E01, ...]] pairs from the calculated energies in this database
        """

        yield from self.get_energies(molecule_name, "%", "%", "%", "%", optimized)

    def get_energies(self, molecule_name, method, basis, cp, tag, optimized = False):
        """
        Returns a generator of pairs of [molecule, energies] where energies is an array of the form [E0, ...] of molecules in the database with the given
        method, basis, cp and tag.

        Includes any optimized geometries.

        % can be used as a wildcard to stand in for any method, basis, cp, or tag.

        Args:
            method  - retrieve only energies computed with this method
            basis   - retrieve only energies computed with this basis
            cp      - retrieve only energies computed with this coutnerpoise correction
            tag     - retrieve only energies with this tag
            optimized - if True, then only retrieve those energies that use an optimized geometry

        Returns:
            a generator of [molecule, [E0, E1, E2, E01, ...]] pairs from the calculated energies in this database using the given model and tag
        """
        
        # get a list of all molecules with the given name
        molecule_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Molecules WHERE name=?", (molecule_name,)).fetchall()]

        # get the id of the model corresponding to the given method, basis, and cp
        model_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Models WHERE method LIKE ? AND basis LIKE ? AND cp LIKE ?", (method, basis, cp)).fetchall()]

        # get a list of all calculations that have the appropriate method, basis, and cp
        calculation_ids = []

        # loop over all the selected molecules
        for molecule_id in molecule_ids:

            # loop over all selected models
            for model_id in model_ids:

                # if optimized is true, get only those calculations which are marked as optimized in the database
                if optimized:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND tag LIKE ? AND optimized=?", (model_id, molecule_id, tag, 1)).fetchall()]

                # otherwise, get all energies (even those marked as optimized)
                else:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND tag LIKE ?", (model_id, molecule_id, tag)).fetchall()]
        
        # loop over all the selected calculations
        for calculation_id in calculation_ids:
            
            # get the molecule id corresponding to this calculation
            molecule_id = self.cursor.execute("SELECT molecule_id FROM Calculations WHERE ROWID=?", (calculation_id,)).fetchone()[0]

            # check if any energies are uncomputed
            is_complete = True
            for job_id in [job_id_tuple[0] for job_id_tuple in self.cursor.execute("SELECT job_id FROM Energies WHERE calculation_id=?", (calculation_id,)).fetchall()]:
                if self.cursor.execute("SELECT status FROM Jobs WHERE ROWID=?", (job_id,)).fetchone()[0] != "completed":
                    is_complete = False
                    break

            # else statement is only ran if inner for loop
            if not is_complete:
                continue

            # get the energies corresponding to this calculation
            energies = [energy_index_value_pair[1] for energy_index_value_pair in sorted(self.cursor.execute("SELECT energy_index, energy FROM Energies WHERE calculation_id=?", (calculation_id,)).fetchall())]

            # Reconstruct the molecule from the information in this database
            molecule = self.get_molecule(molecule_id)
            
            yield molecule, energies

    def count_calculations(self, molecule_name = "%", method = "%", basis = "%", cp = "%", tag = "%", optimized = False):
        """
        Counts the number of calculations in the database with the given molecule, model, and tag

        % can be used as a wildcard to stand in for any method, basis, cp, or tag. 

        Args:
            molecule_name    - count only calculations with this moleucle
            method      - the model's method
            basis       - the model's basis
            cp          - the model's cp
            tag         - count only calculations with this cp
            optimized   - if True, then only count calculations for optimized geometries

        Returns:
            the number of calculations in the database as a 4 item tuple: [pending calculations, partly calculations, completed calculations, failed calculations]
        """

        # get a list of all molecules with the given name
        molecule_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Molecules WHERE name LIKE ?", (molecule_name,)).fetchall()]

        # get the id of the model corresponding to the given method, basis, and cp
        model_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Models WHERE method LIKE ? AND basis LIKE ? AND cp LIKE ?", (method, basis, cp)).fetchall()]

        # get a list of all calculations that have the appropriate method, basis, and cp
        calculation_ids = []

        # loop over all selected molecules
        for molecule_id in molecule_ids:

            # loop over all selected models
            for model_id in model_ids:

                # if optimized is True, only get calculations that use optimized geometries
                if optimized:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND tag LIKE ? AND optimized=?", (model_id, molecule_id, tag, 1)).fetchall()]

                # if optimized is False, get all calculations (even those with optimized geometries)
                else:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND tag LIKE ?", (model_id, molecule_id, tag)).fetchall()]

        # count the number of calculations with each status
        pending = 0
        partly = 0
        completed = 0
        failed = 0

        # loop over each selected calculation
        for calculation_id in calculation_ids:

            # get a list of all the jobs for that calculation
            job_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT job_id FROM Energies WHERE calculation_id=?", (calculation_id,)).fetchall()]

            # count the number of jobs with each status
            num_pending = 0
            num_running = 0
            num_completed = 0
            num_failed = 0

            # loop over all the jobs for this calculation
            for job_id in job_ids:

                # get the status of this job
                status = self.cursor.execute("SELECT status FROM Jobs WHERE ROWID=?", (job_id,)).fetchone()[0]

                # increment the corresponding counter variable
                if status == "pending":
                    num_pending += 1
                elif status == "running":
                    num_running += 1
                elif status == "completed":
                    num_completed += 1
                elif status == "failed":
                    num_failed += 1

            # if any jobs from this calculation failed, the calculation is considered to have failed
            if num_failed > 0:
                failed += 1

            # if no jobs failed, are pending, or are running, then the calculation is complete
            elif num_pending == 0 and num_running == 0:
                completed += 1

            # if no jobs failed, are running, or are completed, then the calculation is pending
            elif num_running == 0 and num_completed == 0:
                pending += 1

            # otherwise the calculation is considered partly done
            else:
                partly += 1

        return pending, partly, completed, failed

    def count_energies(self, molecule_name = "%", method = "%", basis = "%", cp = "%", tag = "%", optimized = False):
        """
        Counts the number of energies in the database with the given molcule, model, and tag

        % can be used as a wildcard to stand in for any method, basis, cp, or tag. 

        Args:
            molecule_name    - count only energies with this moleucle
            method      - the model's method
            basis       - the model's basis
            cp          - the model's cp
            tag         - count only energies with this cp
            optimized   - if True, then only count energies for optimized geometries

        Returns:
            the number of energies in the database as a 4 item tuple: [pending, running, completed, failed]
        """

        # get a list of all molecules with the given name
        molecule_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Molecules WHERE name LIKE ?", (molecule_name,)).fetchall()]

        # get the id of the model corresponding to the given method, basis, and cp
        model_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Models WHERE method LIKE ? AND basis LIKE ? AND cp LIKE ?", (method, basis, cp)).fetchall()]

        # get a list of all calculations that have the appropriate method, basis, and cp
        calculation_ids = []

        # loop over all selected molecules
        for molecule_id in molecule_ids:

            # loop over all selected models
            for model_id in model_ids:

                # if optimized is True, only get calculations that use optimized geometries
                if optimized:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND tag LIKE ? AND optimized=?", (model_id, molecule_id, tag, 1)).fetchall()]

                # if optimized is False, get all calculations (even those with optimized geometries)
                else:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND tag LIKE ?", (model_id, molecule_id, tag)).fetchall()]

        # count the number of energies with each status
        pending = 0
        running = 0
        completed = 0
        failed = 0

        # loop over all selected calculations
        for calculation_id in calculation_ids:

            # get a list of all the jobs for this calculation
            job_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT job_id FROM Energies WHERE calculation_id=?", (calculation_id,)).fetchall()]

            # loop over all the jobs for this calculation
            for job_id in job_ids:

                # get the status of this job
                status = self.cursor.execute("SELECT status FROM Jobs WHERE ROWID=?", (job_id,)).fetchone()[0]

                # increment the corresponding counter variable
                if status == "pending":
                    pending += 1
                elif status == "running":
                    running += 1
                elif status == "completed":
                    completed += 1
                elif status == "failed":
                    failed += 1

        return pending, running, completed, failed

    def what_models(self, molecule_name = "%", tag = "%", optimized = False):
        """
        Given a molecule name, yields information about the models that exist for that molecule and how many calculations are associated with each of them

        Args:
            molecule_name - the name of the molecule to get information about
            tag         - only look at calculations with this tag, '%' serves as a wild-card
            optimized   - if True, only consider calculations which are marked as optimized. Default is False

        Returns:
            A generator object which yields tuples of the form (method, basis, cp, (pending, partly, completed, failed))
        """

        # get a list of all the molecules in the database with the given name
        molecule_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT ROWID FROM Molecules WHERE name LIKE ?", (molecule_name,)).fetchall()]

        model_ids = []
        for molecule_id in molecule_ids:
            model_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT model_id FROM Calculations WHERE molecule_id=? AND tag LIKE ? AND optimized=?", (molecule_id, tag, optimized)).fetchall()]

        model_ids = set(model_ids)

        for model_id in model_ids:
            method, basis, cp = self.cursor.execute("SELECT method, basis, cp FROM Models WHERE ROWID=?", (model_id,)).fetchone()
            calculation_count = self.count_calculations(molecule_name, method, basis, cp, tag, optimized)

            yield method, basis, True if cp == 1 else False, calculation_count

    def what_molecules(self, method = "%", basis = "%", cp = "%", tag = "%", optimized = False):
        """
        Given a model, yields information about the molecules that exist for that model and how many calculations are associated with each of them

        Args:
            method  - the model's method
            basis   - the model's basis
            cp      - the model's cp
            tag     - only consier calculations with this tag
            optimized - if True, only consider calculations marked as optimized geometries

        Returns:
            A generator object which yields tuples of the form (molecule_name, (pending, partly, completed, failed))
        """
        try:
            model_id = self.cursor.execute("SELECT ROWID FROM Models WHERE method LIKE ? AND basis LIKE ? AND cp LIKE ?", (method, basis, cp)).fetchone()[0]
        except TypeError:
            return

        molecule_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute("SELECT molecule_id FROM Calculations WHERE model_id=? AND tag LIKE ? AND optimized=?", (model_id, tag, optimized)).fetchall()]

        molecule_names = []

        for molecule_id in molecule_ids:
            molecule_names.append(self.cursor.execute("SELECT name FROM Molecules WHERE ROWID=?", (molecule_id,)).fetchone()[0])

        molecule_names = set(molecule_names)

        for molecule_name in molecule_names:
            calculation_count = self.count_calculations(molecule_name, method, basis, cp, tag, optimized)

            yield molecule_name, calculation_count
            
    def get_symmetry(self, molecule_name):
        molecule_id = self.cursor.execute("SELECT ROWID FROM Molecules WHERE name=?", (molecule_name,)).fetchone()[0]

        return self.get_molecule(molecule_id).get_symmetry()

    def get_molecule(self, molecule_id):
        # Reconstruct the molecule from the information in the database
        molecule = Molecule()
        
        # loop over all rows in the Fragments table that correspond to this molecule
        for fragment_id, name, charge, spin in self.cursor.execute("SELECT ROWID, name, charge, spin FROM Fragments WHERE molecule_id=?", (molecule_id,)).fetchall():
            fragment = Fragment(name, charge, spin)

            # loop over all rows in the Atoms table that correspond to this fragment
            for symbol, symmetry_class, x, y, z in self.cursor.execute("SELECT symbol, symmetry_class, x, y, z FROM Atoms WHERE fragment_id=?", (fragment_id,)).fetchall():
                fragment.add_atom(Atom(symbol, symmetry_class, x, y, z))

            molecule.add_fragment(fragment)

        return molecule
    

    def clean(self):
        """
        Goes thru all jobs in the database, and sets any that are "running" to "pending". Should only be used when you are prepared to give up on any currently running jobs
        """

        self.cursor.execute("UPDATE Jobs SET status=? WHERE status=?", ("pending", "running"))

class Job(object):
    """
    Contains all the information the user needs to be able to make a calculation
    """
    def __init__(self, molecule, method, basis, cp, fragments, job_id):
        """
        Creates a new Job with the given arguments

        Args:
            molecule - this calculation's molecule
            method  - method to use to run the calculation
            basis   - basis to use to run the calculation
            cp      - whether to use cointerpoise correction for this calculation
            fragments - fragments to include in this calculation, specified as an array of indicies
            job_id  - id of the job corresponding to this energy
        """
        self.molecule = molecule
        self.method = method
        self.basis = basis
        self.cp = cp
        self.fragments = fragments
        self.job_id = job_id

class Calculation(object):
    """
    Contains all the information pertaining to a single calculation
    """
    def __init__(self, molecule, method, basis, cp, tag, energies):
        self.molecule = molecule
        self.method = method
        self.basis = basis
        self.cp = cp
        self.tag = tag
        self.energies = energies

def number_of_energies(number_of_fragments, cp):
    """
    Returns the number of nmer energies that a molecule with number_of_fragments fragments will have
    1 -> 1
    2 -> 3
    3 -> 7
    4 -> 15

    Equal to the size of the power set of a set of size number_of_fragments minus the empty set

    cp corrected models will cause 1 extra energy for each monomer in the molecule

    Args:
        number_of_fragments - the number of fragments in a molecule
        cp                  - if cp is enabled, then 1 energy will be added for each monomer for the non-cp corrected versions of the monomer energies used to computed the deformation energies

    Returns:
        the number of energies a molecule with the given number of fragments will have
    """

    return 2 ** number_of_fragments - 1  + (number_of_fragments if cp and number_of_fragments > 1 else 0)

# THESE TWO METHODS COULD PROBABLY BE BETTER, BUT THEY WORK

def energy_index_to_fragment_indicies(energy_index, number_of_fragments, cp):
    """
    Returns an array of fragment indicies that should be included in a calculation with energy_index index in a molecule with number_of_fragments fragments

    Args:
        energy_index - index of this energy
        number_of_fragments - number of fragments in the molecule that this energy belongs to
        cp - if cp-correction is on. If true then there is one additional energy per monomer, for the non-cp enabled calculations.

    Returns:
        array of the indicies of the fragments included in this energy
    """

    if energy_index < 0 or energy_index >= number_of_energies(number_of_fragments, cp):
        raise ValueError("energy_index", energy_index, "in range [0, {}) for {} fragments with cp={}".format(number_of_energies(number_of_fragments, cp), number_of_fragments, cp))

    combinations = [inner for outer in (itertools.combinations(range(number_of_fragments), combination_size) for combination_size in range(1, number_of_fragments + 1)) for inner in outer]

    if energy_index >= len(combinations) and cp and number_of_fragments > 1:
        return combinations[energy_index - len(combinations)]

    return combinations[energy_index]
