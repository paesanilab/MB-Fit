# external package imports
import itertools, time, datetime, psycopg2, numpy as np, pandas as pd

# absolute module imports
from potential_fitting.utils import files
from potential_fitting.molecule import Atom, Fragment, Molecule
from potential_fitting.exceptions import InconsistentDatabaseError, InvalidValueError


class Database():
    """
    Database class. Allows one to access a database and perform operations on it.
    """
    
    def __init__(self):
        """
        Initializes the database, sets up connection, and cursor.

        Args:

        Returns:
            A new Database object.
        """

        # connection is used to get the cursor, commit to the database, and close the database
        self.connection = psycopg2.connect("host='piggy.pl.ucsd.edu' port=5432 dbname='potential_fitting' user='potential_fitting' password='9t8ARDuN2Wy49VtMOrcJyHtOzyKhkiId'")
        # the cursor is used to execute operations on the database
        self.cursor = self.connection.cursor()

    # the __enter__() and __exit__() methods define a database as a context manager, meaning you can use
    # with ... as ... syntax on it

    def __enter__(self):
        """
        Simply returns self.

        Args:
            None.

        Returms:
            self.
        """

        return self

    def __exit__(self, type, value, traceback):
        """
        Automatically saves and closes the database.

        Args:
            None.

        Returns:
            None.
        """

        self.save()
        self.close()
        
        # returning false lets the context manager know that no exceptions were handled in the __exit__() method
        return False

    def save(self):
        """
        Save any Changes to the database.
        Changes will not be perminent untill this function is called.

        Args:
            None.
        
        Returns:
            None.
        """

        self.connection.commit()

    def close(self):
        """
        Close the database.

        Always close the database after you are done using it.

        Closing the database does NOT automatically save it.

        Args:
            None.

        Returns:
            None.
        """

        self.connection.close()

    def create(self):
        """
        Creates all the tables in the database, if they do not exist already.

        Args:
            None.

        Returns:
            None.
        """

    def exists(self, table, **property_value_pairs):
        if len(property_value_pairs) < 1:
            # must specify at least one property
            raise Exception
        properties = tuple(property_value_pairs)
        properties_string = properties[0] + "=%s"
        for property in properties[1:]:
            properties_string += " AND "
            properties_string += property + "=%s"

        values = tuple(property_value_pairs.values())

        self.cursor.execute("SELECT EXISTS(SELECT * FROM {} WHERE {})".format(table, properties_string), values)

        return self.cursor.fetchone()[0]

    def insert(self, table, **property_value_pairs):

        properties = tuple(property_value_pairs)

        property_string = "(" + properties[0]
        for prop in properties[1:]:
            property_string += ", " + prop
        property_string += ")"

        values = tuple(property_value_pairs.values())

        values_string = "(%s"
        for value in values[1:]:
            values_string += ", %s"
        values_string += ")"

        self.cursor.execute("INSERT INTO {}{} VALUES {}".format(table, property_string, values_string), values)
        return None

    def select(self, table, fetch_all, *values, **property_value_pairs):
        if len(property_value_pairs) < 1:
            # must specify at least one property
            raise Exception
        if len(values) < 1:
            # must specify at least one value
            raise Exception

        properties = tuple(property_value_pairs)
        properties_string = properties[0] + "=%s"
        for property in properties[1:]:
            properties_string += " AND "
            properties_string += property + "=%s"

        query_vals = tuple(property_value_pairs.values())

        if values == ("*",):
            values_string = "*"
            values = ()
        else:
            values_string = values[0]
            for value in values[1:]:
                values_string += ", " + value
            values_string += ""

        self.cursor.execute("SELECT {} FROM {} WHERE {}".format(values_string, table, properties_string), query_vals)

        if fetch_all:
            return self.cursor.fetchall()

        return self.cursor.fetchone()

    def update(self, table, values, **property_value_pairs):
        pass

    def create_postgres_array(self, *values):
        return "{" + ",".join([str(i) for i in values]) + "}"

    def create_current_postgres_time(self):
        offset = time.timezone if (time.localtime().tm_isdst == 0) else time.altzone
        offset = offset / 60 / 60 * -1

        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + str(int(offset))

    def add_atom_info(self, atom):
        """
        Adds a single atom type's info to the atom_info table if it does not already exist
        """
        if not self.exists("atom_info", atomic_symbol = atom.get_name()):
            self.insert("atom_info", atomic_symbol = atom.get_name())

    def add_fragment_info(self, fragment):

        # use transactions

        """
        Adds a single fragment type's info to the fragment_info table if it does not already exist
        """
        if not self.exists("fragment_info", name = fragment.get_name(), charge = fragment.get_charge(), spin = fragment.get_spin_multiplicity()):
            # updates fragment_info table
            self.insert("fragment_info", name = fragment.get_name(), charge = fragment.get_charge(), spin = fragment.get_spin_multiplicity())

            # updates atom_info table
            for atom in fragment.get_atoms():
                self.add_atom_info(atom)

            # updates fragment_contents
            atoms = [[atom.get_name(), atom.get_symmetry_class()] for atom in fragment.get_atoms()]

            symbol_symmetry_pairs, counts = np.unique(atoms, return_counts = True, axis = 0)
            counts = [int(i) for i in counts]

            for symbol_symmetry_pair, count in zip(symbol_symmetry_pairs, counts):
                atomic_symbol = symbol_symmetry_pair[0]
                symmetry_symbol = symbol_symmetry_pair[1]
                self.insert("fragment_contents", frag_name = fragment.get_name(), atom_symbol = atomic_symbol, count = count, symmetry = symmetry_symbol)

    def add_molecule_info(self, molecule):

        """
        Adds a single molecule type's info to the molecule_info table if it does not already exist
        """

        if not self.exists("molecule_info", name = molecule.get_name()):
            self.insert("molecule_info", name = molecule.get_name())

            for fragment in molecule.get_fragments():
                self.add_fragment_info(fragment)

            fragments, counts = np.unique([fragment.get_name() for fragment in molecule.get_fragments()], return_counts = True)
            counts = [int(i) for i in counts]

            for fragment, count in zip(fragments, counts):
                self.insert("molecule_contents", mol_name = molecule.get_name(), frag_name = fragment, count = count)

    def add_model_info(self, method, basis, cp):
        model ="{}/{}/{}".format(method, basis, cp)
        if not self.exists("model_info", name = model):
            self.insert("model_info", name = model)

    def add_molecule(self, molecule):
        # should always be at COM and Principle axes before calling!!!!

        if not self.exists("molecule_list", mol_hash = molecule.get_SHA1()):

            coordinates = []

            for fragment in molecule.get_fragments():
                for atom in fragment.get_atoms():
                    coordinates.append(atom.get_x())
                    coordinates.append(atom.get_y())
                    coordinates.append(atom.get_z())

            self.insert("molecule_list", molecule.get_SHA1(), molecule.get_name(), coordinates)

            self.add_molecule_info(molecule)

    def add_calculation(self, molecule, method, basis, cp, *tags):
        model ="{}/{}/{}".format(method, basis, cp)
        if not self.exists("molecule_properties", mol_hash = molecule.get_SHA1(), model_name = model):

            self.add_molecule(molecule)
            self.add_model_info(method, basis, cp)

            for n in range(1, molecule.get_num_fragments() + 1):
                for fragment_indices in itertools.combinations(range(molecule.get_num_fragments()), n):
                    self.insert("molecule_properties", molecule.get_SHA1(), model, list(fragment_indices), [], [], "pending", [], None)
            
            if not self.exists("tags", mol_hash = molecule.get_SHA1(), model_name = model):
                self.insert("tags", molecule.get_SHA1(), model, [])


        for tag in tags:
            self.cursor.execute("SELECT EXISTS(SELECT * FROM tags WHERE mol_hash=%s AND model_name=%s AND %s=ANY(tag_names))", (molecule.get_SHA1(), model, tag))
            if not self.cursor.fetchone()[0]:

                self.cursor.execute("UPDATE tags SET tag_names=(%s || tag_names) WHERE mol_hash=%s AND model_name=%s", (self.create_postgres_array(tag), molecule.get_SHA1(), model))

    def get_molecule(self, mol_hash):
        mol_name, atom_coordinates = self.select("molecule_list", False, "mol_name", "atom_coordinates", mol_hash = mol_hash)

        molecule = Molecule()

        molecule_contents = self.select("molecule_contents", True, "frag_name", "count", mol_name = mol_name)

        for frag_name, frag_count in molecule_contents:
            for i in range(frag_count):
                charge, spin = self.select("fragment_info", False, "charge", "spin", name = frag_name)

                fragment = Fragment(frag_name, charge, spin)

                fragment_contents = self.select("fragment_contents", True, "atom_symbol", "symmetry", "count", frag_name = frag_name)
                for atom_symbol, atom_symmetry, atom_count in fragment_contents:
                    for k in range(atom_count):
                        atom = Atom(atom_symbol, atom_symmetry, 0, 0, 0)

                        fragment.add_atom(atom)

                molecule.add_fragment(fragment)

        for atom in molecule.get_atoms():
            atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
            atom_coordinates = atom_coordinates[3:]

        return molecule

    def get_calculation(self, client_name, *tags):
        try:
            mol_hash, model_name, frag_indices = self.select("molecule_properties", False, "mol_hash", "model_name", "frag_indices", status = "pending")
        except TypeError:
            raise Exception

        self.insert("log_files", start_time = self.create_current_postgres_time() ,client_name = client_name)

        self.cursor.execute("SELECT last_value FROM log_files_log_id_seq")
        log_id = self.cursor.fetchone()[0]

        self.cursor.execute("UPDATE molecule_properties SET status=%s, most_recent_log_id=%s WHERE mol_hash=%s AND model_name=%s AND frag_indices=%s", ("dispatched", log_id, mol_hash, model_name, self.create_postgres_array(*frag_indices)))
        molecule = self.get_molecule(mol_hash)


        return molecule, model_name, frag_indices

    def set_properties(self, molecule, model_name, frag_indices, energy, log_file):
        # possibly check to make sure molecule is not morphed

        self.cursor.execute("UPDATE molecule_properties SET energies=%s, status=%s WHERE mol_hash=%s AND model_name=%s AND frag_indices=%s", (self.create_postgres_array(energy), "complete", molecule.get_SHA1(), model_name, self.create_postgres_array(*frag_indices)))
        log_id = self.select("molecule_properties", False, "most_recent_log_id", mol_hash = molecule.get_SHA1(), model_name = model_name, frag_indices = self.create_postgres_array(*frag_indices))

        # make sure to push previous log id to past_log_ids

        self.cursor.execute("UPDATE log_files SET end_time=%s, log_text=%s WHERE log_id=%s", (self.create_current_postgres_time(), log_file, log_id))

    def reset_all_calculations(self, *tag):
        self.cursor.execute("UPDATE molecule_properties SET status=%s WHERE status=%s", ("pending", "dispatched"))
        #temporary
        self.cursor.execute("UPDATE molecule_properties SET status=%s WHERE status=%s", ("pending", "complete"))

    #------------------------------------------------------------------------------------#
    #------------------------------------------OLD DATABASE CODE BELOW THIS LINE---------#
    #------------------------------------------------------------------------------------#
    '''

    def add_calculation(self, molecule, method, basis, cp, *tags, optimized = False):
        """
        Adds a calculation to the database. A calculation consists of pending nmer energies to be calculated in a given
        method, basis, and cp.

        If cp is False, then all monomer, dimer, trimer, etc energies will be added to the database. If cp is True,
        then monomer energies will be added in their own basis set in addition to their counterpoise corrected basis
        sets.

        Args:
            molecule        - The molecule associated with this calculation.
            method          - The method to use to calculate the energies.
            basis           - The basis to use to calculate the energies.
            cp              - Whether to use counterpoise correction for this calculation.
            tags            - A set of tags to be able to identify this calculation at a later date.
            optimized       - A special flag which indicates whether this calculation uses an optimized geometry.

        Returns:
            None.
        """

        # check if this model is not already in Models table
        if not self.cursor.execute("SELECT EXISTS(SELECT * FROM Models WHERE method=? AND basis=? AND cp=?)", (method,
                basis, cp)).fetchone()[0]:

            # create entry in Models table
            self.cursor.execute("INSERT INTO Models (method, basis, cp) VALUES (?, ?, ?)", (method, basis, cp))
            
            # get id of this model
            model_id = self.cursor.lastrowid

        else:

            # get id of this model
            model_id = self.cursor.execute("SELECT ROWID FROM Models WHERE method=? AND basis=? AND cp=?", (method,
                    basis, cp)).fetchone()[0]

        # get the SHA1 hash of this molecule
        molecule_hash = molecule.get_SHA1()

        # check if this molecule is not already in Molecules table
        if not self.cursor.execute("SELECT EXISTS(SELECT * FROM Molecules WHERE name=? AND hash=?)",
                (molecule.get_name(), molecule_hash)).fetchone()[0]:

            # create entry in Molecules table
            self.cursor.execute("INSERT INTO Molecules (name, hash) VALUES (?, ?)", (molecule.get_name(),
                    molecule_hash))
        
            # get id of this molecule
            molecule_id = self.cursor.lastrowid

            # insert molecule's fragments into the table
            for fragment in molecule.get_fragments():
                self.cursor.execute("INSERT INTO Fragments (molecule_id, frag_index, name, charge, spin) VALUES (?, ?, ?, ?, ?)",
                        (molecule_id, fragment.get_index(), fragment.get_name(), fragment.get_charge(), fragment.get_spin_multiplicity()))
                
                # get id of this fragment
                fragment_id = self.cursor.lastrowid

                # insert fragment's atoms into the table
                for atom in fragment.get_atoms():
                   self.cursor.execute(
                            "INSERT INTO Atoms (fragment_id, symbol, symmetry_class, x, y, z) VALUES "
                            + "(?, ?, ?, ?, ?, ?)", (fragment_id, atom.get_name(), atom.get_symmetry_class(),
                            atom.get_x(), atom.get_y(), atom.get_z())) 

        else:
            
            # get id of this molecule
            molecule_id = self.cursor.execute("SELECT ROWID FROM Molecules WHERE name=? AND hash=?",
                    (molecule.get_name(), molecule_hash)).fetchone()[0]

        # check if the calculation is not already in the Calculations table

        # alphabetize the tags
        tag_string = " ".join(sorted(set(tags)))

        if not self.cursor.execute(
                "SELECT EXISTS(SELECT * FROM Calculations WHERE molecule_id=? AND model_id=? AND tag LIKE ? AND "
                + "optimized=?)", (molecule_id, model_id, "%", optimized)).fetchone()[0]:
            
            # create entry in Calculations table
            self.cursor.execute("INSERT INTO Calculations (molecule_id, model_id, tag, optimized) VALUES (?, ?, ?, ?)",
                    (molecule_id, model_id, tag_string, optimized))

            # get id of this calculation
            calculation_id = self.cursor.lastrowid

            # add rows to the Energies table
            for energy_index in range(number_of_energies(molecule.get_num_fragments(), cp)):
                # create a job for this energy
                self.cursor.execute("INSERT INTO Jobs (status) VALUES (?)", ("pending",))

                # get the id of this job
                job_id = self.cursor.lastrowid
        
                # insert row into Energies table for this energy
                self.cursor.execute("INSERT INTO Energies (calculation_id, energy_index, job_id) VALUES (?, ?, ?)",
                        (calculation_id, energy_index, job_id))

        else:
            calc_id, existing_tags = self.cursor.execute(
                    "SELECT ROWID, TAG FROM Calculations WHERE molecule_id=? AND model_id=? AND tag LIKE ? AND "
                    + "optimized=?", (molecule_id, model_id, "%", optimized)).fetchone()

            tags = tags + tuple(existing_tags.split(" "))
            tag_string = " ".join(sorted(set(tags)))

            self.cursor.execute("UPDATE Calculations SET tag=? WHERE ROWID=?", (tag_string, calc_id))

    def get_missing_energy(self):
        """
        Returns a single Job object, which contains the info a user needs to calculate a pending energy in the table.

        The user should calculate the energy, then call either set_energy() if the calculation converged, or 
        set_failed() if it did not.

        The Job will be updating from pending to running.

        Args:
            None.

        Returns:
            A Job object containing all the info associated with a calculation to be performed. Will instead return
            None if there are no more energies to be calculated in the database.
        """

        # retrieve a pending job from the Jobs table
        try:
            job_id = self.cursor.execute("SELECT ROWID FROM Jobs WHERE status=?", ("pending",)).fetchone()[0]
        except TypeError:
            # this error will be thrown if there are no more energies to calculate, because None[0] throws a TypeError
            return None
        
        # retrieve the calculation and energy to be calculated for this job from the Energies table
        calculation_id, energy_index = self.cursor.execute(
                "SELECT calculation_id, energy_index FROM Energies WHERE job_id=?", (job_id,)).fetchone()

        # retrieve the molecule and model from the Calculations table
        molecule_id, model_id, tags, optimized = self.cursor.execute(
                "SELECT molecule_id, model_id, tag, optimized FROM Calculations WHERE ROWID=?",
                (calculation_id,)).fetchone()

        # retrieve the method, basis, and cp for this model
        method, basis, cp = self.cursor.execute("SELECT method, basis, cp FROM Models WHERE ROWID=?",
                (model_id,)).fetchone()

        # Reconstruct the molecule from the information in this database
        molecule = self.get_molecule(molecule_id)

        # get the indices of the fragments to include in this calculation
        fragment_indices = energy_index_to_fragment_indices(energy_index, molecule.get_num_fragments(),
                True if cp == 1 else False)

        # update the job with its start date and running status
        self.cursor.execute("UPDATE Jobs SET status=?, start_date=? WHERE ROWID=?", ("running",
                datetime.datetime.today().strftime('%Y/%m/%d'), job_id))
        
        return Job(molecule, method, basis, True if cp == 1 and energy_index <
                number_of_energies(molecule.get_num_fragments(), False) else False, fragment_indices, job_id)

    def missing_energies(self):
        """
        Generates Jobs for all the pending energies in the database. Each Job contains all the information needed
        to perform one calculation.

        The user should calculate each energy, then call either set_energy() if the calculation converged, or 
        set_failed() if it did not.

        Args:
            None.

        Returns:
            Yields Job objects until one has been generated for every pending job in the database.
        """
        while True:
            calculation = self.get_missing_energy()
            
            if calculation is None:
                return

            yield calculation

    def set_energy(self, job_id, energy, log_file):
        """
        Updates the database by filling in the energy for one Job.

        The Job will be updated from running to completed.

        Args:
            job_id          - The id of the Job associated with this energy, included in the Job object retrieved by
                    get_missing_energy() or missing_energies().
            energy          - The calculated energy for the given Job in Hartrees.
            log_file        - Local path to the log file for the given Job.

        Returns:
            None.
        """
        
        if self.cursor.execute("SELECT status FROM Jobs WHERE ROWID=?", (job_id,)).fetchone()[0] != "running":
            print("This job is not set to running, but tried to have its energy set. Most likely, this job was "
                    + "dispatched a long time ago, and the database has since given up on ever getting a response. "
                    + "Energy will be placed into the database anyways.")


        # update the row in the Energies table corresponding to this energy
        self.cursor.execute("UPDATE Energies SET energy=? WHERE job_id=?", (energy, job_id))

        # update the information about this job
        self.cursor.execute("UPDATE Jobs SET status=?, log_file=?, end_date=? WHERE ROWID=?", ("completed", log_file,
                datetime.datetime.today().strftime('%Y/%m/%d'), job_id))

    def set_failed(self, job_id, log_path):
        """
        Notifies the database that a specific Job failed to have its energy calculated. The energy could have not
        converged, or failed for any number of reasons.

        Args:
            job_id          - The id of the job which failed, included in the Job object retrieved by
                    get_missing_energy() or missing_energies().
            log_file        - Local path to the log file for the given Job.

        Returns:
            None.
        """

        self.cursor.execute("UPDATE Jobs SET status=?, log_file=? WHERE ROWID=?", ("failed", log_path, job_id))

    def get_complete_energies(self, molecule_name, optimized = False):
        """
        Generates pairs of [molecule, energies] where energies is an array of the form [E0, E1, E2, E01, ...].

        If optimized is true, this will only generate pairs corresponding to optimized geometries. Otherwise, it will
        generate all pairs, including those with optimized geometries.

        Args:
            molecule_name   - The name of the molecule to get complete energies for.
            optimized       - If True, only get optimized energies.

        Returns:
            Yields [molecule, [E0, E1, E2, E01, ...]] pairs from the calculated energies in this database.
        """

        yield from self.get_energies(molecule_name, "%", "%", "%", "%", optimized)

    def get_energies(self, molecule_name, method, basis, cp, *tags, optimized = False):
        """
        Generates pairs of [molecule, energies] where energies is an array of the form [E0, E1, E2, E01, ...].

        Only gets those energies computed using the given method, basis, cp, and marked with at least 1 of the given
        tags.

        If optimized is true, this will only generate pairs corresponding to optimized geometries. Otherwise, it will
        generate all pairs, including those with optimized geometries.

        % can be used as a wildcard to stand in for any method, basis, cp, or tag.

        Args:
            method          - Retrieve only energies computed with this method.
            basis           - Retrieve only energies computed with this basis.
            cp              - Retrieve only energies computed with this coutnerpoise correction.
            tags            - Retrieve only energies marked with at least one of the given tags.
            optimized       - If True, then only retrieve those energies that are associated with an optimized
                    geometry.

        Returns:
            Yields [molecule, [E0, E1, E2, E01, ...]] pairs from the calculated energies in this database using the
            given method, basis, cp, and marked with at least one of the given tags.
        """
        
        # get a list of all molecules with the given name
        molecule_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                "SELECT ROWID FROM Molecules WHERE name=?", (molecule_name,)).fetchall()]

        # get the id of the model corresponding to the given method, basis, and cp
        model_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                "SELECT ROWID FROM Models WHERE method LIKE ? AND basis LIKE ? AND cp LIKE ?",
                (method, basis, cp)).fetchall()]

        # get a list of all calculations that have the appropriate method, basis, and cp
        calculation_ids = []

        # loop over all the selected molecules
        for molecule_id in molecule_ids:

            # loop over all selected models
            for model_id in model_ids:

                if len(tags) == 0:
                    tag_like_operation = "TRUE"
                else:
                    tag_like_operation = "(" + " OR ".join(["tag LIKE '{0} %' OR tag LIKE '% {0} %' OR tag LIKE '% {0}' OR tag='{0}'".format(tag) for tag in tags]) + ")"
                print(tag_like_operation)

                # if optimized is true, get only those calculations which are marked as optimized in the database
                if optimized:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                            "SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND " 
                            + tag_like_operation + " AND optimized=?", (model_id, molecule_id, 1)).fetchall()]

                # otherwise, get all energies (even those marked as optimized)
                else:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                            "SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND "
                            + tag_like_operation, (model_id, molecule_id)).fetchall()]
        
        # loop over all the selected calculations
        for calculation_id in calculation_ids:
            
            # get the molecule id corresponding to this calculation
            molecule_id = self.cursor.execute("SELECT molecule_id FROM Calculations WHERE ROWID=?",
                    (calculation_id,)).fetchone()[0]

            # check if any energies are uncomputed
            is_complete = True
            for job_id in [job_id_tuple[0] for job_id_tuple in self.cursor.execute(
                    "SELECT job_id FROM Energies WHERE calculation_id=?", (calculation_id,)).fetchall()]:
                if self.cursor.execute("SELECT status FROM Jobs WHERE ROWID=?",
                        (job_id,)).fetchone()[0] != "completed":
                    is_complete = False
                    break

            # else statement is only ran if inner for loop
            if not is_complete:
                continue

            # get the energies corresponding to this calculation
            energies = [energy_index_value_pair[1] for energy_index_value_pair in sorted(self.cursor.execute(
                    "SELECT energy_index, energy FROM Energies WHERE calculation_id=?", (calculation_id,)).fetchall())]

            # Reconstruct the molecule from the information in this database
            molecule = self.get_molecule(molecule_id)
            
            yield molecule, energies

    def count_calculations(self, molecule_name = "%", method = "%", basis = "%", cp = "%", tags = ["%"],
            optimized = False):
        """
        Counts the number of calculations in the database with the given molecule, method, basis, and at least one of
        the given tags.

        % can be used as a wildcard to stand in for any molecule name, method, basis, cp, or tag. 

        Args:
            molecule_name   - Count only calculations with a molecule by this name.
            method          - Count only calculations with this method.
            basis           - Count only calculations with this basis.
            cp              - Count only calculations iwth this cp.
            tags            - Count only calculations marked with at least one of these tags.
            optimized       - If True, then only count calculations for optimized geometries. Otherwise, counts
                    calculations for both optimized and unoptimized geometries.

        Returns:
            The number of calculations in the database as a 4 item tuple: [pending calculations, partly calculations,
            completed calculations, failed calculations] A partly calculation is a calculation whose energies are a
            mixture of pending, running, completed, and failed, or all running.
        """

        # get a list of all molecules with the given name
        molecule_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                "SELECT ROWID FROM Molecules WHERE name LIKE ?", (molecule_name,)).fetchall()]

        # get the id of the model corresponding to the given method, basis, and cp
        model_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                "SELECT ROWID FROM Models WHERE method LIKE ? AND basis LIKE ? AND cp LIKE ?", (method, basis,
                cp)).fetchall()]

        # get a list of all calculations that have the appropriate method, basis, and cp
        calculation_ids = []

        # loop over all selected molecules
        for molecule_id in molecule_ids:

            # loop over all selected models
            for model_id in model_ids:

                if len(tags) == 0:
                    tag_like_operation = "TRUE"
                else:
                    tag_like_operation = "(" + " OR ".join(["tag LIKE '{0} %' OR tag LIKE '% {0} %' OR tag LIKE '% {0}' OR tag='{0}'".format(tag) for tag in tags]) + ")"

                # if optimized is true, get only those calculations which are marked as optimized in the database
                if optimized:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                            "SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND " 
                            + tag_like_operation + " AND optimized=?", (model_id, molecule_id, 1)).fetchall()]

                # otherwise, get all energies (even those marked as optimized)
                else:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                            "SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND "
                            + tag_like_operation, (model_id, molecule_id)).fetchall()]
        
        # count the number of calculations with each status
        pending = 0
        partly = 0
        completed = 0
        failed = 0

        # loop over each selected calculation
        for calculation_id in calculation_ids:

            # get a list of all the jobs for that calculation
            job_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                    "SELECT job_id FROM Energies WHERE calculation_id=?", (calculation_id,)).fetchall()]

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

    def count_energies(self, molecule_name = "%", method = "%", basis = "%", cp = "%", tags = ["%"],
            optimized = False):
        """
        Counts the number of energies in the database with the given molecule, method, basis, and at least one of the
        give tags.

        % can be used as a wildcard to stand in for any molecule name, method, basis, cp, or tag.

        Differes from count_calculations(), because count_energies() will count each nmer energy of a calculation,
        while count_calculations() will only count once for each calculation.

        Args:
            molecule_name   - Count only energies with a molecule by this name.
            method          - Count only eneriges with this method.
            basis           - Count only energies with this basis.
            cp              - Count only energies with this cp.
            tags            - Count only calculations marked with at least one of these tags.
            optimized       - If True, then only count energies for optimized geometries. Otherwise, counts
                    energies for both optimized and unoptimized geometries.

        Returns:
            The number of energies in the database as a 4 item tuple: [pending energies, running energies,
            completed energies, failed energies]. 
        """

        # get a list of all molecules with the given name
        molecule_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                "SELECT ROWID FROM Molecules WHERE name LIKE ?", (molecule_name,)).fetchall()]

        # get the id of the model corresponding to the given method, basis, and cp
        model_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                "SELECT ROWID FROM Models WHERE method LIKE ? AND basis LIKE ? AND cp LIKE ?", (method, basis,
                cp)).fetchall()]

        # get a list of all calculations that have the appropriate method, basis, and cp
        calculation_ids = []

        # loop over all selected molecules
        for molecule_id in molecule_ids:

            # loop over all selected models
            for model_id in model_ids:

                if len(tags) == 0:
                    tag_like_operation = "TRUE"
                else:
                    tag_like_operation = "(" + " OR ".join(["tag LIKE '{0} %' OR tag LIKE '% {0} %' OR tag LIKE '% {0}' OR tag='{0}'".format(tag) for tag in tags]) + ")"

                # if optimized is true, get only those calculations which are marked as optimized in the database
                if optimized:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                            "SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND " 
                            + tag_like_operation + " AND optimized=?", (model_id, molecule_id, 1)).fetchall()]

                # otherwise, get all energies (even those marked as optimized)
                else:
                    calculation_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                            "SELECT ROWID FROM Calculations WHERE model_id=? AND molecule_id=? AND "
                            + tag_like_operation, (model_id, molecule_id)).fetchall()]
        
        # count the number of energies with each status
        pending = 0
        running = 0
        completed = 0
        failed = 0

        # loop over all selected calculations
        for calculation_id in calculation_ids:

            # get a list of all the jobs for this calculation
            job_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                    "SELECT job_id FROM Energies WHERE calculation_id=?", (calculation_id,)).fetchall()]

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

    def what_models(self, molecule_name = "%", tags = ["%"], optimized = False):
        """
        Given a molecule name, yields information about the models that exist for that molecule and how many
        calculations are associated with each of them.

        Args:
            molecule_name   - The name of the molecule to get information about.
            tags            - Only look at calculations with at least one of these tags, '%' serves as a wild-card
            optimized       - If True, only consider calculations which are marked as optimized. Otherwise, consider
                    both optimized and unoptimized calculations.

        Returns:
            Yields tuples of the form (method, basis, cp, (pending, partly, completed,
            failed)), where pending, partly, completed, failed correspond to the number of calculations with each
            status for the corresponding method, basis, cp for the queried molecule.
        """

        # get a list of all the molecules in the database with the given name
        molecule_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                "SELECT ROWID FROM Molecules WHERE name LIKE ?", (molecule_name,)).fetchall()]

        model_ids = []
        for molecule_id in molecule_ids:

            if len(tags) == 0:
                tag_like_operation = "TRUE"
            else:
                tag_like_operation = "(" + " OR ".join(["tag LIKE '{0} %' OR tag LIKE '% {0} %' OR tag LIKE '% {0}' OR tag='{0}'".format(tag) for tag in tags]) + ")"

            model_ids += [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                    "SELECT model_id FROM Calculations WHERE molecule_id=? AND "+ tag_like_operation+" AND optimized=?",
                    (molecule_id, optimized)).fetchall()]

        model_ids = set(model_ids)

        for model_id in model_ids:
            method, basis, cp = self.cursor.execute("SELECT method, basis, cp FROM Models WHERE ROWID=?",
                    (model_id,)).fetchone()
            calculation_count = self.count_calculations(molecule_name, method, basis, cp, tags, optimized)

            yield method, basis, True if cp == 1 else False, calculation_count

    def what_molecules(self, method = "%", basis = "%", cp = "%", tags = ["%"], optimized = False):
        """
        Given a model, yields information about the molecules that exist for that model and how many calculations are 
        associated with each of them.

        Args:
            method          - The model's method.
            basis           - The model's basis.
            cp              - The model's cp.
            tags            - Only consier calculations with at least one of these tags.
            optimized       - If True, only consider calculations marked as optimized geometries. Otherwise, consider
                    both optimized and unoptimized calculations.

        Returns:
            Yields tuples of the form (molecule_name, (pending, partly, completed,
            failed)), where pending, partly, completed, failed correspond to the number of calculations with each
            status for the corresponding molecule_name.
        """
        try:
            model_id = self.cursor.execute(
                    "SELECT ROWID FROM Models WHERE method LIKE ? AND basis LIKE ? AND cp LIKE ?", (method, basis,
                    cp)).fetchone()[0]
        except TypeError:
            return

        if len(tags) == 0:
            tag_like_operation = "TRUE"
        else:
            tag_like_operation = "(" + " OR ".join(["tag LIKE '{0} %' OR tag LIKE '% {0} %' OR tag LIKE '% {0}' OR tag='{0}'".format(tag) for tag in tags]) + ")"

        molecule_ids = [fetch_tuple[0] for fetch_tuple in self.cursor.execute(
                "SELECT molecule_id FROM Calculations WHERE model_id=? AND " + tag_like_operation + " AND optimized=?",
                (model_id, optimized)).fetchall()]

        molecule_names = []

        for molecule_id in molecule_ids:
            molecule_names.append(self.cursor.execute("SELECT name FROM Molecules WHERE ROWID=?",
                    (molecule_id,)).fetchone()[0])

        molecule_names = set(molecule_names)

        for molecule_name in molecule_names:
            calculation_count = self.count_calculations(molecule_name, method, basis, cp, tags, optimized)

            yield molecule_name, calculation_count
            
    def get_symmetry(self, molecule_name):
        molecule_id = self.cursor.execute("SELECT ROWID FROM Molecules WHERE name=?", (molecule_name,)).fetchone()[0]

        return self.get_molecule(molecule_id).get_symmetry()

    def get_molecule(self, molecule_id):
        # Reconstruct the molecule from the information in the database
        molecule = Molecule()

        index_fragment_pairs = []
        
        # loop over all rows in the Fragments table that correspond to this molecule
        for fragment_id, index, name, charge, spin in self.cursor.execute(
                "SELECT ROWID, frag_index, name, charge, spin FROM Fragments WHERE molecule_id=?", (molecule_id,)).fetchall():
            fragment = Fragment(name, charge, spin)

            # loop over all rows in the Atoms table that correspond to this fragment
            for symbol, symmetry_class, x, y, z in self.cursor.execute(
                    "SELECT symbol, symmetry_class, x, y, z FROM Atoms WHERE fragment_id=?",
                    (fragment_id,)).fetchall():
                fragment.add_atom(Atom(symbol, symmetry_class, x, y, z))

            index_fragment_pairs.append([index, fragment])

        for index, fragment in sorted(index_fragment_pairs, key=lambda pair: pair[0]):
            molecule.add_fragment(fragment)

        return molecule
    

    def clean(self):
        """
        Goes through all jobs in the database, and sets any that are "running" to "pending". Should only be used when 
        you are prepared to give up on any currently running jobs.

        Args:
            None.

        Returns:
            None.
        """

        self.cursor.execute("UPDATE Jobs SET status=? WHERE status=?", ("pending", "running"))
    '''

class Job(object):
    """
    Contains all the information the user needs to be able to calculate a single energy.
    """
    def __init__(self, molecule, method, basis, cp, fragments, job_id):
        """
        Creates a new Job with the given arguments.

        Args:
            molecule        - This Job's molecule.
            method          - Method to use to run the energy calculation.
            basis           - Basis to use to run the energy calculation.
            cp              - Whether to use cointerpoise correction for this energy calculation.
            fragments       - Fragments to include in this energy calculation, specified as an array of indices.
            job_id          - Id of the job corresponding to this energy. Will be needed to call set_energy().

        Returns:
            A new Job object.
        """

        self.molecule = molecule
        self.method = method
        self.basis = basis
        self.cp = cp
        self.fragments = fragments
        self.job_id = job_id

class Calculation(object):
    """
    Contains all the information pertaining to a single completed Calculation.
    """
    def __init__(self, molecule, method, basis, cp, tags, energies):
        """
        Creates a new Calculation with the given arguments.

        Args:
            molecule        - The molecule of this Calculation.
            method          - The method used to calculate the energies in this Calculation.
            basis           - The basis used to calculate the energies in this Calculation.
            cp              - Whether counterpoise correction was used to calculated the energies in this Calculation.
            tags            - This Calculation's tags.
            energies        - The nmer energies of this calculation as a list [E0, E1, E2, E01, ...].

        Returns:
            A new Calculation object.
        """

        self.molecule = molecule
        self.method = method
        self.basis = basis
        self.cp = cp
        self.tags = tags
        self.energies = energies

def number_of_energies(number_of_fragments, cp):
    """
    Returns the number of nmer energies that a molecule with number_of_fragments fragments will have.
    1 -> 1 [E0]
    2 -> 3 [E0, E1, E01]
    3 -> 7 [E0, E1, E2, E01, E02, E12, E012]
    4 -> 15 etc

    counterpoise corrected models will cause 1 extra energy for each monomer in the molecule. Unless the molecule is a
    monomer, in which case counterpoise correction has no affect.

    Args:
        number_of_fragments - The number of fragments in the molecule.
        cp                  - If cp is enabled, then 1 additional energy will be added for each monomer for the 
                non-cp corrected versions of the monomer energies. These will be used to compute deformation energies.

    Returns:
        The number of nmer energies a molecule with the given number of fragments will have.
    """

    return 2 ** number_of_fragments - 1  + (number_of_fragments if cp and number_of_fragments > 1 else 0)

def energy_index_to_fragment_indices(energy_index, number_of_fragments, cp):
    """
    Returns an array of fragment indices that should be included in a energy calculation with energy_index index in a
    molecule with number_of_fragments fragments.

    The order is first all monomers, then dimers, then trimers, etc. Finally, if cp is set to True, the the non-cp
    corrected monomers are last.

    Args:
        energy_index        - The index of this energy.
        number_of_fragments - The number of fragments in the molecule that this energy belongs to.
        cp                  - If cp-correction is on. If true then there is one additional energy per monomer, for the
                non-cp enabled calculations.

    Returns:
        List of the indices of the fragments included in this energy.
    """

    if energy_index < 0 or energy_index >= number_of_energies(number_of_fragments, cp):
        raise ValueError("energy_index", energy_index, "in range [0, {}) for {} fragments with cp={}".format(number_of_energies(number_of_fragments, cp), number_of_fragments, cp))

    combinations = [inner for outer in (itertools.combinations(range(number_of_fragments), combination_size) for combination_size in range(1, number_of_fragments + 1)) for inner in outer]

    if energy_index >= len(combinations) and cp and number_of_fragments > 1:
        return combinations[energy_index - len(combinations)]

    return combinations[energy_index]
