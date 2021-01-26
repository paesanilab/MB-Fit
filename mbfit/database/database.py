# external package imports
import itertools, numpy as np, copy, sys, os

# absolute module imports
from mbfit.molecule import Atom, Fragment, Molecule
from mbfit.exceptions import PotentialFittingError, NoSuchMoleculeError, DatabaseOperationError, \
        DatabaseInitializationError, DatabaseNotEmptyError, DatabaseConnectionError, InvalidValueError, \
        NoPendingCalculationsError, StandardOrderError, LibraryNotAvailableError
from mbfit.utils import SettingsReader

# only import psycopg2 if it is installed.
try:
    import psycopg2
except ModuleNotFoundError:
    pass

class Database():

    """
    Database class. Allows one to access a database and perform operations on it.
    """
    
    def __init__(self, config_file, batch_size=100):
        """
        Initializer for database object. Opens connection and sets up cursor.

        Args:
            config_file     - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
            batch_size      - number of operations to perfrom on the database per round trip to the server.
                    larger numbers will be more efficient, but you should not exceed a couple thousand.
                    Default is 100.

        Returns:
            A new Database object.
        """

        # Check if psycopg2 is installed.
        try:
            import psycopg2
        except ModuleNotFoundError:
            raise LibraryNotAvailableError("psycopg2")

        self.batch_size = 0
        self.set_batch_size(batch_size)

        # parse the user's config file to get their login info

        config = SettingsReader(config_file)

        host = config.get("database", "host")
        port = config.get("database", "port")
        database = config.get("database", "database")
        username = config.get("database", "username")
        password = config.get("database", "password")

        self.name = host + " " + database

        """
        host = "piggy.pl.ucsd.edu"
        port = "5432"
        database = "mbfit"
        username = "mbfit"
        password = "9t8ARDuN2Wy49VtMOrcJyHtOzyKhkiId"
        """

        # connection is used to get the cursor, commit to the database, and close the database
        try:
            self.connection = psycopg2.connect("host='{}' port={} dbname='{}' user='{}' password='{}'".format(host, port, database, username, password))
        except psycopg2.OperationalError as e:
            raise DatabaseConnectionError(self.name, str(e))

        # the cursor is used to execute operations on the database
        self.cursor = self.connection.cursor()

    # the __enter__() and __exit__() methods define a database as a context manager, meaning you can use
    # with ... as ... syntax on it

    def __enter__(self):
        """
        Simply returns self. Called when entering the context manager.
        Args:
            None.
        Returns:
            self.
        """

        return self

    def __exit__(self, exception, value, traceback):
        """
        Saves and closes the database. Called when exiting the context manager.
        Args:
            None.
        Returns:
            None.
        """

        # only save if all operations executed without error
        if exception is None:
            self.save()

        self.close()

        # returning false lets the context manager know that no exceptions were handled in the __exit__() method
        return False

    def save(self):
        """
        Commits changes to the database.
        Changes will not be permanent until this function is called.
        Args:
            None.

        Returns:
            None.
        """

        self.connection.commit()

    def close(self):
        """
        Closes the database.
        Always close the database after you are done using it.
        Closing the database does NOT automatically save it.
        Any non-saved changes will be lost.
        Args:
            None.
        Returns:
            None.
        """

        self.connection.close()

    def set_batch_size(self, batch_size):
        """
        Sets the number of operations to do on the database per server round trip.
        Larger numbers will be more efficient, but you should not exceed a couple thousand.
        Args:
            batch_size      - The number of operations to do per server round trip.
        Returns:
            None.
        """

        if batch_size < 1:
            raise InvalidValueError("batch_size", batch_size, "must be at least 1.")

        self.batch_size = batch_size

    def get_batch_size(self):
        """
        Gets the batch size, the number of operations performed per server round trip.
        Args:
            None.
        Returns:
            batch_size
        """

        return self.batch_size

    def get_notices(self):
        """
        Gets a list of all notices received from the Database. A notice is a logging message or warning that does
        not crash the program.
        Args:
            None.
        Returns:
            A list of all notices that have been received since this Database object was created.
        """

        return self.connection.notices

    def get_last_notice(self):
        """
        Gets the last notice received from the Database. A notice is a logging message or warning that does not crash
        the program.
        Args:
            None.
        Returns:
            The last notice received from the Database.
        """

        if len(self.connection.notices) > 0:
            return self.connection.notices[-1]
        else:
            return None

    def clear_notices(self):
        """
        Clears the list of notices for this Database. A notice is a logging message or warning that does not crash
        the program.
        Args:
            None.
        Returns:
            None.
        """
        self.connection.notices = []

    def create(self):
        """
        Creates all the tables, procedures, and sequences in the database.
        Args:
            None.
        Returns:
            None.
        """

        self.single_execute("SELECT * FROM pg_catalog.pg_tables WHERE schemaname != %s AND schemaname != %s;",
                            ("pg_catalog", "information_schema"))

        table_names = [table[1] for table in self.cursor.fetchall()]

        if len(table_names) > 0:
            raise DatabaseNotEmptyError(self.name, table_names)

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "init.sql")) as sql_initilizer:
            sql_script = sql_initilizer.read()
            try:
                self.single_execute(sql_script, ())
            except psycopg2.OperationalError as e:
                raise DatabaseInitializationError(self.name, str(e))

    def annihilate(self, confirm="no way"):
        """
        DELETES ALL CONTENT IN ALL TABLES IN THE DATABASE.
        Don't EVER do this for databases you share with others unless you ask them first.
        TODO: Make this method require administrator priveleges on the database.
        In order to protect users, confirm must be specified as
        "confirm" or the content will not be deleted.
        Args:
            confirm         - set to "confirm" if deletion of all content in the database is desired.
        Returns:
            None
        """

        if confirm == "confirm":
            self.execute("PERFORM annihilate();", ())
        else:
            print(
                "annihilate failed. specify confirm = \"confirm\" if deletion of all content in the database is desired.")

    def execute(self, command, params):
        """
        Executes a PostgreSQL command in the database.
        The String command can contain '%s' substrings, each %s
        will be substituted for an item in params. Order of the '%s's
        matches the order in which the params will be substituted.
        One param per '%s'.
        pyscopg2 is pretty good about converting python types to
        the correct postgres types, except lists. For lists, add
        create_postgres_array(list) to params.
        This method is used for executing long commands that must be wrapped in a
        function block.
        Args:
            command         - The command to run.
            params          - Parameters for the command.
        Returns:
            None.
        """

        self.single_execute(
            "DO $$" +
            "   BEGIN " +
            command +
            "END $$",
            params)

    def single_execute(self, command, params):
        """
        Executes a PostgreSQL command in the database.
        The String command can contain '%s' substrings, each %s
        will be substituted for an item in params. Order of the '%s's
        matches the order in which the params will be substituted.
        One param per '%s'.
        pyscopg2 is pretty good about converting python types to
        the correct postgres types, except lists. For lists, add
        create_postgres_array(list) to params.
        This method is used for executing single select commands, insert commands,
        and other simple commands.
        Args:
            command         - The command to run.
            params          - Parameters for the command.
        Returns:
            None.
        """

        try:
            print(command)
            self.cursor.execute(
                command,
                params)
        except psycopg2.OperationalError as e:
            raise DatabaseOperationError(self.name, str(e.diag.message_primary)) from None
        except psycopg2.InternalError as e:
            self.connection.rollback()
            raise DatabaseOperationError(self.name, str(e.diag.message_primary)) from None
        except psycopg2.ProgrammingError as e:
            self.connection.rollback()
            raise DatabaseOperationError(self.name, str(e.diag.message_primary)) from None

    def create_postgres_array(self, *values):
        """
        Creates a postgres array with an arbitrary number of values.
        You should always call this before passing an array as a param
        into execute().
        Args:
            values          - Set of values which will form the elements
                    of the postgres array.
        Returns:
            The postgres array as a python string, ready to be passed into
            execute() as a param.
        """

        return "{" + ",".join([str(i) for i in values]) + "}"

    def add_calculations(self, molecule_list, method, basis, cp, *tags, optimized=False):
        """
        Adds new calculations to the database.
        Will queue the database to calculate the energies of each molecule in the list
        with the given model. Does not calculate the energies.
        All molecules must be of the same type.
        Args:
            molecule_list   - List of molecules whose energies are wanted.
            method          - Method to use to calculate the molecules' energies.
            basis           - Basis to use to calculate the molecules' energies.
            cp              - True if counterpoise correction should be used in the calculation of the molecules' energies.
            tags            - Set of tags to label these calculations in the database.
            optimized       - True if all molecules represent optimized geometries.
        Returns:
            None.
        """

        command_string = ""
        params = []

        batch_count = 0

        order, frag_order, SMILES = None, None, None

        for molecule in molecule_list:
            if order is None:
                order, frag_order = molecule.get_standard_order_order()
                SMILES = [frag.get_standard_SMILE() for frag in molecule.get_standard_order()]

            molecule = molecule.get_reordered_copy(order, frag_order, SMILES)
            coordinates = []
            for fragment in molecule.get_fragments():
                for atom in fragment.get_atoms():
                    coordinates.append(atom.get_x())
                    coordinates.append(atom.get_y())
                    coordinates.append(atom.get_z())

            command_string += "PERFORM add_calculation(%s, %s, ARRAY["
            params += (molecule.get_SHA1(), molecule.get_name())

            fragments = [fragment.get_name() for fragment in molecule.get_fragments()]

            frag_names, counts = np.unique(fragments, return_counts=True, axis=0)
            fragment_counts = [int(i) for i in counts]

            for frag_name in frag_names:

                fragment = None

                for frag in molecule.get_fragments():
                    if frag.get_name() == frag_name:
                        fragment = frag

                atoms = [[atom.get_name(), atom.get_symmetry_class()] for atom in fragment.get_atoms()]

                symbol_symmetry_pairs, counts = np.unique(atoms, return_counts=True, axis=0)
                atom_counts = [int(i) for i in counts]

                symbol_symmetry_count_pairs = [[symbol_symmetry_pairs[i][0], symbol_symmetry_pairs[i][1], counts[i]]
                                               for i in range(len(atom_counts))]
                symbol_symmetry_count_pairs.sort(key=lambda x: x[1])

                symbols = [symbol for symbol, symmetry, count in symbol_symmetry_count_pairs]
                symmetries = [chr(65 + index) for index in range(len(symbol_symmetry_count_pairs))]
                counts = [count for symbol, symmetry, count in symbol_symmetry_count_pairs]
                command_string += "construct_fragment(%s, %s, %s, %s, %s, %s, %s)"
                if not frag_name == frag_names[-1]:
                    command_string += ", "
                params += (frag_name, fragment.get_charge(), fragment.get_spin_multiplicity(), fragment.get_SMILE(),
                           self.create_postgres_array(*symbols),
                           self.create_postgres_array(*symmetries), self.create_postgres_array(*counts))

            command_string += "], %s, %s, %s, %s, %s, %s, %s);"
            params += (
            fragment_counts, coordinates, method, basis, cp, self.create_postgres_array(*tags), optimized)

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    def build_empty_molecule(self, mol_name):
        """
        Returns a copy of the mol_name molecule from inside the database with all atom
        coordinates set to 0.
        Params:
            mol_name        - The name of the molecule to construct.
        Returns:
            a Molecule object with all coordinates set to 0.
        """

        fragments = []

        next_start_symmetry = 'A'

        self.single_execute("SELECT * FROM get_empty_molecule(%s)", (mol_name,))

        frag_info = sorted(self.cursor.fetchall(), key = lambda x: x[0])

        for frag_name, charge, spin, atomic_symbols, atomic_symmetries, atomic_counts, SMILE, frag_count in frag_info:

            for i in range(frag_count):

                next_symmetry = next_start_symmetry

                atoms = []

                fragment_contents = [(atom, symmetry, count) for atom, symmetry, count in
                                     zip(atomic_symbols, atomic_symmetries, atomic_counts)]

                fragment_contents.sort(key=lambda x: x[1])

                for atom_symbol, atom_symmetry, atom_count in fragment_contents:
                    for k in range(atom_count):
                        atom = Atom(atom_symbol, next_symmetry, 0, 0, 0)

                        atoms.append(atom)

                    next_symmetry = chr(ord(next_symmetry) + 1)

                atoms.sort(key=lambda x: x.get_symmetry_class())

                fragments.append(Fragment(atoms, frag_name, charge, spin, SMILE))

            next_start_symmetry = next_symmetry

        molecule = Molecule(fragments)

        for index, atom in enumerate(molecule.get_atoms()):
            atom.set_xyz(index, index, index)

        return molecule

    def count_pending_calculations(self, *tags):
        self.single_execute("SELECT * FROM count_pending_calculations(%s)", (self.create_postgres_array(*tags),))

        count = self.cursor.fetchone()[0]

        return count

    def count_dispatched_calculations(self):
        self.single_execute("SELECT * FROM count_dispatched_calculations()", ())

        count = self.cursor.fetchone()[0]

        return count

    def get_all_calculations(self, client_name, *tags, calculations_to_do=sys.maxsize):
        """
        Gets uncalculaed energies from the database so that the user can calculate them.
        Pass the output into set_properties to update the energies in the database.
        Args:
            client_name     - The name of the client that will perform these calculations.
            tags            - Only fetch calculations with these tags.
            calculations_to_do - Maximum number of calculations to fetch. Defualt is unlimited.
        Yields:
            (molecule, method, basis, cp, use_cp, frag_indices)
            molecule        - The molecule whose energy should be calculated.
            method          - Method to do the calculation.
            basis           - Basis to do the calculation.
            cp              - True if the model uses counterpoise correction.
            use_cp          - True if counterpoise correction should be used for this calculation.
            frag_indices    - List of indices of fragments that should be included in the calculation.
                    If use_cp is True, then include other fragments as ghost atoms.
            cp is not the same as use_cp. Some models have cp, but should not
            use cp for some of their energies.
        """

        name_to_order_dict = {}

        while True:

            self.single_execute("SELECT * FROM get_pending_molecule_name(%s)", (self.create_postgres_array(*tags),))

            molecule_name = self.cursor.fetchone()[0]

            if molecule_name == "":
                break

            self.single_execute("SELECT * FROM get_pending_calculations(%s, %s, %s, %s)", (
            molecule_name, client_name, self.create_postgres_array(*tags),
            min(self.batch_size, calculations_to_do)))

            pending_calcs = self.cursor.fetchall()

            calculations_to_do -= len(pending_calcs)

            empty_molecule = self.build_empty_molecule(molecule_name)
            empty_molecule = empty_molecule.get_standard_copy()

            for atom_coordinates, model, frag_indices, use_cp in pending_calcs:
                molecule = copy.deepcopy(empty_molecule)

                for atom in molecule.get_atoms():
                    atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
                    atom_coordinates = atom_coordinates[3:]

                method = model[:model.index("/")]
                model = model[model.index("/") + 1:]
                basis = model[:model.index("/")]
                cp = model[model.index("/") + 1:] == "True"

                try:
                    order, frag_orders, SMILES = name_to_order_dict[molecule.get_name()]
                except KeyError:
                    order, frag_orders = molecule.get_standard_order_order()
                    SMILES = [frag.get_standard_SMILE() for frag in molecule.get_standard_order()]
                    name_to_order_dict[molecule.get_name()] = order, frag_orders, SMILES

                molecule = molecule.get_reordered_copy(order, frag_orders, SMILES)

                yield molecule, method, basis, cp, use_cp, frag_indices

            if calculations_to_do < 1:
                return

    def set_properties(self, calculation_results, overwrite=False):
        """
        Sets newly calculated energies in the database.
        Args:
            calculation_results - List of tuples of format
                    (molecule, method, basis, cp, use_cp, frag_indices, result, energy, log_text)
                    molecule    - Molecule whose energy should be set.
                    method      - Method used to calculate the new energy.
                    basis       - Basis used to calculate the new energy.
                    cp          - True if the model for this energy includes counterpoise correction.
                    use_cp      - True if counterpoise correction was used for this calculation.
                    frag_indices - Fragments included in this calculation.
                    result      - True if the calculation succeeded.
                    energy      - The calculated energy in atomic units. Not used in result is False.
                    log_text    - The text of the log file for this calculation.
        Returns:
            None.
        """

        command_string = ""
        params = []

        batch_count = 0

        name_to_order_dict = {}

        for molecule, method, basis, cp, use_cp, frag_indices, result, energy, log_text in calculation_results:
            try:
                order, frag_orders, SMILES = name_to_order_dict[molecule.get_name()]
            except KeyError:
                order, frag_orders = molecule.get_standard_order_order()
                SMILES = [frag.get_standard_SMILE() for frag in molecule.get_standard_order()]
                name_to_order_dict[molecule.get_name()] = order, frag_orders, SMILES

            molecule = molecule.get_reordered_copy(order, frag_orders, SMILES)

            model_name = method + "/" + basis + "/" + str(cp)

            command_string += "PERFORM set_properties(%s, %s, %s, %s, %s, %s, %s, %s);"
            params += [molecule.get_SHA1(), model_name, use_cp, self.create_postgres_array(*frag_indices), result,
                       energy, log_text, overwrite]

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    def get_1B_training_set(self, molecule_name, names, SMILES, method, basis, cp, *tags):
        """
        Gets a 1B training set from the calculated energies in the database.
        All complete calculations which match the given method, basis, cp, and tags
        will be included.
        Args:
            molecule_name   - Name of the molecule for which a training set is desired.
            names           - List of name of the monomer.
            SMILES          - List of SMILE string of the monomer, the atoms in the training set will
                    be in this order.
            method          - Method of this training set.
            basis           - Basis of this training set.
            cp              - Counterpoise correction of this training set.
            tags            - Only include calculations marked with at least one of these tags.
        Yields:
            (molecule, energy)
            molecule        - One molecule in the training set.
            energy          - Its 1B deformation energy.
        """

        model_name = "{}/{}/{}".format(method, basis, cp)
        batch_offset = 0

        self.single_execute("SELECT * FROM count_entries(%s)", (molecule_name,))
        max_count = self.cursor.fetchone()[0]

        empty_molecule = self.build_empty_molecule(molecule_name)

        order, frag_orders = None, None

        while True:
            self.single_execute("SELECT * FROM get_1B_training_set(%s, %s, %s, %s, %s)", (
            molecule_name, model_name, self.create_postgres_array(*tags), batch_offset, self.batch_size))
            training_set = self.cursor.fetchall()

            for atom_coordinates, energy in training_set:
                molecule = copy.deepcopy(empty_molecule)

                for atom in molecule.get_atoms():
                    atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
                    atom_coordinates = atom_coordinates[3:]

                if order is None:
                    order, frag_orders = molecule.get_reorder_order(names, SMILES)

                yield molecule.get_reordered_copy(order, frag_orders, SMILES), energy

            batch_offset += self.batch_size

            if batch_offset > max_count:
                return

    def get_training_set_size(self, names, method, basis, cp, *tags):
        model_name = "{}/{}/{}".format(method, basis, cp)

        standard_names = sorted(names)

        molecule_name = "-".join(standard_names)

        self.single_execute("SELECT * FROM count_training_set_size(%s, %s, %s)",
                            (molecule_name, model_name, self.create_postgres_array(*tags)))

        count = self.cursor.fetchone()[0]

        return count

    def get_training_set(self, names, SMILES, method, basis, cp, *tags):

        self.clear_notices()

        model_name = "{}/{}/{}".format(method, basis, cp)
        batch_offset = 0

        order, frag_orders, energies_order = None, None, None

        standard_names = sorted(names)

        molecule_name = "-".join(standard_names)

        max_count = self.get_training_set_size(names, method, basis, cp, *tags)

        empty_molecule = self.build_empty_molecule(molecule_name)

        while True:
            self.single_execute("SELECT * FROM get_training_set(%s, %s, %s, %s, %s, %s)", (
                molecule_name, self.create_postgres_array(*standard_names), model_name,
                self.create_postgres_array(*tags), batch_offset,
                self.batch_size))
            training_set = self.cursor.fetchall()

            for atom_coordinates, binding_energy, nb_energy, deformation_energies in training_set:
                molecule = copy.deepcopy(empty_molecule)

                for atom in molecule.get_atoms():
                    atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
                    atom_coordinates = atom_coordinates[3:]

                if order is None:
                    order, frag_orders = molecule.get_reorder_order(names, SMILES)
                    energies_order = Database.get_energies_order(order, molecule.get_num_fragments(), False)

                deformation_energies = [deformation_energies[i] for i in energies_order[:len(deformation_energies)]]

                yield molecule.get_reordered_copy(order, frag_orders,
                                                  SMILES), binding_energy, nb_energy, deformation_energies

            batch_offset += self.batch_size

            if batch_offset > max_count:                
                if self.get_last_notice() is not None and "Multiple optimized geometries" in self.get_last_notice():
                    print(self.get_last_notice(), "Using the lowest energy optimized geometry to calculate deformation"
                                                              " energies for this training set.")
                return

    def get_2B_training_set(self, molecule_name, names, SMILES, method, basis, cp, *tags):
        """
        Gets a 2B training set from the calculated energies in the database.
        All complete calculations which match the given method, basis, cp, and tags
        will be included.
        Args:
            molecule_name   - Name of the molecule for which a training set is desired.
            names           - List of names of the two monomers, the training set will have the monomers
                    in this order.
            SMILES          - List of SMILE strings of each monomer, the atoms in the training set will
                    be in this order.
            method          - Method of this training set.
            basis           - Basis of this training set.
            cp              - Counterpoise correction of this training set.
            tags            - Only include calculations marked with at least one of these tags.
        Yields:
            (molecule, binding_energy, interaction_energy, monomer1_energy, monomer2_energy)
            molecule        - One molecule in the training set.
            binding_energy  - Its 2B binding energy.
            interaction_energy - Its 2B interaction energy.
            monomer1_energy - Its first monomer's deformation energy.
            monomer2_energy - Its second monomer's deformation energy.
        """

        model_name = "{}/{}/{}".format(method, basis, cp)
        batch_offset = 0

        order, frag_orders = None, None

        monomer1_name, monomer2_name = sorted([names[0], names[1]])

        molecule_name = monomer1_name + "-" + monomer2_name

        self.single_execute("SELECT * FROM count_entries(%s)", (molecule_name,))
        max_count = self.cursor.fetchone()[0]

        empty_molecule = self.build_empty_molecule(molecule_name)

        while True:
            self.single_execute("SELECT * FROM get_2B_training_set(%s, %s, %s, %s, %s, %s, %s)", (
            molecule_name, monomer1_name, monomer2_name, model_name, self.create_postgres_array(*tags),
            batch_offset, self.batch_size))
            training_set = self.cursor.fetchall()

            for atom_coordinates, binding_energy, interaction_energy, monomer1_energy, monomer2_energy in training_set:
                molecule = copy.deepcopy(empty_molecule)

                for atom in molecule.get_atoms():
                    atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
                    atom_coordinates = atom_coordinates[3:]

                if order is None:
                    order, frag_orders = molecule.get_reorder_order(names, SMILES)

                if order == [1, 0]:
                    monomer1_energy, monomer2_energy = monomer2_energy, monomer1_energy

                yield molecule.get_reordered_copy(order, frag_orders,
                                                  SMILES), binding_energy, interaction_energy, monomer1_energy, monomer2_energy

            batch_offset += self.batch_size

            if batch_offset > max_count:
                return

    def export_calculations(self, names, SMILES, method, basis, cp, *tags):

        model_name = "{}/{}/{}".format(method, basis, cp)
        batch_offset = 0

        order, frag_orders, energies_order = None, None, None

        molecule_name = "-".join(sorted(names))

        self.single_execute("SELECT * FROM count_entries(%s)", (molecule_name,))
        max_count = self.cursor.fetchone()[0]

        empty_molecule = self.build_empty_molecule(molecule_name)

        while True:
            self.single_execute("SELECT * FROM export_calculations(%s, %s, %s, %s, %s)", (
                molecule_name, model_name, self.create_postgres_array(*tags), batch_offset, self.batch_size))
            calculations = self.cursor.fetchall()

            for atom_coordinates, energies in calculations:
                molecule = copy.deepcopy(empty_molecule)

                for atom in molecule.get_atoms():
                    atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
                    atom_coordinates = atom_coordinates[3:]

                if order is None:
                    order, frag_orders = molecule.get_reorder_order(names, SMILES)
                    energies_order = self.get_energies_order(order, molecule.get_num_fragments(), cp)

                yield molecule.get_reordered_copy(order, frag_orders,
                                                  SMILES), [energies[i] for i in energies_order]

            batch_offset += self.batch_size

            if batch_offset > max_count:
                return

    def import_calculations(self, molecule_energies_pairs, method, basis, cp, *tags, optimized=False):
        """
        Imports already completed calculations into the database.
        Args:
            molecule_energies_pairs - Array of 2-tuples (molecule, nmer_energies)
                    each pair being a molecule and its energies in the numeric order
                    by fragment index. If cp is true, cp energies should come before non-cp energies.
            method              - Method used to calculate these energies.
            basis               - Basis used to calculate these energies.
            cp                  - True if counterpoise correction was used for these energies.
            tags                - Tags to mark these energies with.
            optimized           - True if all geometries in molecule_energy_pairs are optimized in
                    the given method and basis.
        """

        command_string = ""
        params = []

        batch_count = 0
        order, frag_order, SMILES, energies_order = None, None, None, None

        for molecule, energies in molecule_energies_pairs:
            if order is None:
                order, frag_order = molecule.get_standard_order_order()
                SMILES = [frag.get_standard_SMILE() for frag in molecule.get_standard_order()]
                energies_order = Database.get_energies_order(order, molecule.get_num_fragments(), cp)

            molecule = molecule.get_reordered_copy(order, frag_order, SMILES)

            energies = [energies[index] for index in energies_order]

            coordinates = []
            for fragment in molecule.get_fragments():
                for atom in fragment.get_atoms():
                    coordinates.append(atom.get_x())
                    coordinates.append(atom.get_y())
                    coordinates.append(atom.get_z())

            command_string += "PERFORM import_calculation(%s, %s, ARRAY["
            params += (molecule.get_SHA1(), molecule.get_name())

            fragments = [fragment.get_name() for fragment in molecule.get_fragments()]

            frag_names, counts = np.unique(fragments, return_counts=True, axis=0)
            fragment_counts = [int(i) for i in counts]

            for frag_name in frag_names:

                fragment = None

                for frag in molecule.get_fragments():
                    if frag.get_name() == frag_name:
                        fragment = frag

                atoms = [[atom.get_name(), atom.get_symmetry_class()] for atom in fragment.get_atoms()]

                symbol_symmetry_pairs, counts = np.unique(atoms, return_counts=True, axis=0)
                atom_counts = [int(i) for i in counts]

                symbols = [symbol for symbol, symmetry in symbol_symmetry_pairs]
                symmetries = [symmetry for symbol, symmetry in symbol_symmetry_pairs]

                command_string += "construct_fragment(%s, %s, %s, %s, %s, %s, %s)"
                if not frag_name == frag_names[-1]:
                    command_string += ", "
                params += (frag_name, fragment.get_charge(), fragment.get_spin_multiplicity(), fragment.get_SMILE(),
                           self.create_postgres_array(*symbols),
                           self.create_postgres_array(*symmetries), self.create_postgres_array(*counts))

            command_string += "], %s, %s, %s, %s, %s, %s, %s, %s);"
            params += (
            fragment_counts, coordinates, method, basis, cp, self.create_postgres_array(*tags), optimized,
            self.create_postgres_array(*energies))

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    def get_failed(self, molecule_name, names, SMILES, method, basis, cp, *tags, optimized=False):
        """
        Gets geometries of energy caclulations that have failed.
        All failed calculations which match the given method, basis, cp, and tags
        will be included.
        Args:
            molecule_name   - Name of the molecule to find failed calculations for.
            names           - List of names of the two monomers, the molecules will have the monomers
                    in this order.
            SMILES          - List of SMILE strings of each monomer, the atoms in the fragments will
                    be in this order.
            method          - Method for the failed calculations.
            basis           - Basis for the failed calculations.
            cp              - Counterpoise correction for the failed calculations.
            tags            - Only include calculations marked with at least one of these tags.
        Yields:
            (molecule, energy, used_cp)
            molecule        - One molecule in the training set.
            frag_indices    - The indices of the fragments used in the failed calculation.
            used_cp         - True of fragments not in frag_indices where included as ghost
                    atoms in this failed calculation.
        """

        model_name = "{}/{}/{}".format(method, basis, cp)
        batch_offset = 0

        self.single_execute("SELECT * FROM count_entries(%s)", (molecule_name,))
        max_count = self.cursor.fetchone()[0]

        empty_molecule = self.build_empty_molecule(molecule_name)

        order, frag_orders = None, None

        while True:
            self.single_execute("SELECT * FROM get_failed_configs(%s, %s, %s, %s, %s)", (
                molecule_name, model_name, self.create_postgres_array(*tags), batch_offset, self.batch_size))
            training_set = self.cursor.fetchall()

            for atom_coordinates, frag_indices, used_cp in training_set:
                molecule = copy.deepcopy(empty_molecule)

                for atom in molecule.get_atoms():
                    atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
                    atom_coordinates = atom_coordinates[3:]

                if order is None:
                    order, frag_orders = molecule.get_reorder_order(names, SMILES)

                yield molecule.get_reordered_copy(order, frag_orders, SMILES), frag_indices, used_cp

            batch_offset += self.batch_size

            if batch_offset > max_count:
                return

    def reset_all_calculations(self, *tags):
        """
        Resets all dispatched, complete, and failed calculations in the database back to pending. Their
        energies are queued for recalculation.
        Args:
            tags            - Currently unused.
        Returns:
            None.
        """
        self.single_execute("UPDATE molecule_properties SET status=%s WHERE status=%s", ("pending", "dispatched"))
        self.single_execute("UPDATE molecule_properties SET status=%s WHERE status=%s", ("pending", "complete"))
        self.single_execute("UPDATE molecule_properties SET status=%s WHERE status=%s", ("pending", "failed"))

    def reset_dispatched(self, *tags):
        """
        Resets all dispatched calculations in the database back to pending. Their
        energies are queued for recalculation.
        Args:
            tags            - Currently unused.
        Returns:
            None.
        """
        self.execute("PERFORM reset_dispatched(%s);", [self.create_postgres_array(*tags)])

    def reset_failed(self, *tags):
        """
        Resets all failed calculations in the database back to pending. Their
        energies are queued for recalculation.
        Args:
            tags            - Currently unused.
        Returns:
            None.
        """
        self.execute("PERFORM reset_failed(%s);", [self.create_postgres_array(*tags)])

    def delete_calculations(self, molecule_list, method, basis, cp, *tags, delete_complete_calculations=False):
        """
        Removes the specified tags from any calculations in the database that matches one of the molecules in
        molecule_list and the given method, basis, and cp.
        Will never delete completed or dispatched energies unless delete_complete_calculations is True, only remove tags from them.
        Will fully delete uncomplete calculations from the database.
        Args:
            molecule_list - Remove tags from calculations involving these molecules.
            method  - Remove tags from calculations with this method.
            basis   - Remove tags from calculations with this basis.
            cp      - Remove tags from calculations with this cp.
            tags    - The tags to remove.
            delete_complete_calculations - If True, delete calculations even if their energy is already
                calculated.
        Returns:
            None.
        """

        command_string = ""
        params = []

        batch_count = 0

        order, frag_order, SMILES = None, None, None

        for molecule in molecule_list:
            if order is None:
                order, frag_order = molecule.get_standard_order_order()
                SMILES = [frag.get_standard_SMILE() for frag in molecule.get_standard_order()]

            molecule = molecule.get_reordered_copy(order, frag_order, SMILES)

            command_string += "PERFORM delete_calculation(%s, %s, %s, %s, %s, %s, %s);"
            params += [molecule.get_SHA1(), molecule.get_name(), method, basis, cp,
                       self.create_postgres_array(*tags), delete_complete_calculations]

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    def delete_all_calculations(self, molecule_name, method, basis, cp, *tags, delete_complete_calculations=False):
        """
        Removes tags from molecules in the database that match the molecule_name.
        Will never delete completed or dispatched energies unless delete_complete_calculations is True, only remove tags from them.
        Will fully delete uncomplete calculations from the database.
        Args:
            molecule_name - Only delete molecules with this name.
            method  - Remove tags from calculations with this method.
            basis   - Remove tags from calculations with this basis.
            cp      - Remove tags from calculations with this cp.
            tags    - The tags to remove.
            delete_complete_calculations - If True, delete calculations even if their energy is already
                calculated.
        """
        self.execute("PERFORM delete_all_calculations(%s, %s, %s, %s, %s, %s);", (
        molecule_name, method, basis, cp, self.create_postgres_array(*tags), delete_complete_calculations))

    def grant_admin_privilege(self, username, *tags):
        """
        Grants admin privileges on a training set to user username.
        Current user must have admin privileges to perform this action.
        Admins can add/remove other admins, writers, and readers.
        Args:
            username    - Username to grant admin privileges to.
            tags        - All training sets to grant admin privileges on.
        Returns:
            None.
        """

        command_string = ""
        params = []

        batch_count = 0

        for tag in tags:

            command_string += "PERFORM grant_admin_privilege(%s, %s);"
            params += (tag, username)

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    def grant_write_privilege(self, username, *tags):
        """
        Grants write privileges on a training set to user username.
        Current user must have admin privileges to perform this action.
        Writers can add/remove elements from the training set.
        Args:
            username    - Username to grant write privileges to.
            tags        - All training sets to grant write privileges on.
        Returns:
            None.
        """

        command_string = ""
        params = []

        batch_count = 0

        for tag in tags:

            command_string += "PERFORM grant_write_privilege(%s, %s);"
            params += (tag, username)

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    def grant_read_privilege(self, username, *tags):
        """
        Grants read privileges on a training set to user username.
        Current user must have admin privileges to perform this action.
        Readers can look at the training set but not modify it.
        Args:
            username    - Username to grant read privileges to.
            tags        - All training sets to grant read privileges on.
        Returns:
            None.
        """

        command_string = ""
        params = []

        batch_count = 0

        for tag in tags:

            command_string += "PERFORM grant_read_privilege(%s, %s);"
            params += (tag, username)

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    def revoke_admin_privilege(self, username, *tags):
        """
        Revokes admin privileges on a training set from user username.
        Current user must have admin privileges to perform this action.
        Admins can add/remove other admins, writers, and readers.
        Args:
            username    - Username to revoke admin privileges from.
            tags        - All training sets to revoke admin privileges from.
        Returns:
            None.
        """

        command_string = ""
        params = []

        batch_count = 0

        for tag in tags:

            command_string += "PERFORM revoke_admin_privilege(%s, %s);"
            params += (tag, username)

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    def revoke_write_privilege(self, username, *tags):
        """
        Revokes write privileges on a training set from user username.
        Current user must have admin privileges to perform this action.
        Writers can add/remove elements from the training set.
        Args:
            username    - Username to revoke write privileges from.
            tags        - All training sets to revoke write privileges from.
        Returns:
            None.
        """

        command_string = ""
        params = []

        batch_count = 0

        for tag in tags:

            command_string += "PERFORM revoke_write_privilege(%s, %s);"
            params += (tag, username)

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    def revoke_read_privilege(self, username, *tags):
        """
        Revokes read privileges on a training set from user username.
        Current user must have admin privileges to perform this action.
        Readers can look at the training set but not modify it.
        Args:
            username    - Username to revoke read privileges from.
            tags        - All training sets to revoke read privileges from.
        Returns:
            None.
        """

        command_string = ""
        params = []

        batch_count = 0

        for tag in tags:

            command_string += "PERFORM revoke_read_privilege(%s, %s);"
            params += (tag, username)

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    @staticmethod
    def get_energies_order(order, num_bodies, cp):
        pre_energy_order = list(range(len(Database.get_permutations(num_bodies, cp))))

        pre_energy_indices = [Database.energy_index_to_frag_indices(i, num_bodies, cp) for i in pre_energy_order]

        post_energy_indices = [(tuple(sorted(order[i] for i in pre)), cp) for pre, cp in pre_energy_indices]

        post_energy_order = [Database.frag_indices_to_energy_index(i, num_bodies, cp) for i in post_energy_indices]

        return post_energy_order

    @staticmethod
    def energy_index_to_frag_indices(index, num_bodies, cp):

        return Database.get_permutations(num_bodies, cp)[index]

    @staticmethod
    def frag_indices_to_energy_index(frag_indices, num_bodies, cp):

        return Database.get_permutations(num_bodies, cp).index(frag_indices)

    @staticmethod
    def get_permutations(num_bodies, cp):
        permutations = []
        for i in range(1, num_bodies + 1):
            perms = list(itertools.combinations(range(num_bodies), i))
            for p in perms:
                if cp and i == 1 and i < num_bodies:
                    permutations.append((p, True))
                permutations.append((p, False))

        return permutations
