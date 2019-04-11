# external package imports
import itertools, psycopg2, numpy as np, copy, sys, os

# absolute module imports
from potential_fitting.molecule import Atom, Fragment, Molecule
from potential_fitting.exceptions import InconsistentDatabaseError, InvalidValueError, NoPendingCalculationsError

class Database():

    """
    Database class. Allows one to access a database and perform operations on it.
    """
    
    def __init__(self, batch_size = 100):
        """
        Initializes the database, sets up connection, and cursor.

        Args:
            batch_size      - number of operations to perfrom on the database per round trip to the server.
                    larger numbers will be more efficient, but you should not exceed a couple thousand.

        Returns:
            A new Database object.
        """

        self.batch_size = batch_size

        # connection is used to get the cursor, commit to the database, and close the database
        #self.connection = psycopg2.connect("host='piggy.pl.ucsd.edu' port=5432 dbname='potential_fitting' user='potential_fitting' password='9t8ARDuN2Wy49VtMOrcJyHtOzyKhkiId'")
        self.connection = psycopg2.connect("host='piggy.pl.ucsd.edu' port=5432 dbname='potential_fitting' user='potential_fitting' password='9t8ARDuN2Wy49VtMOrcJyHtOzyKhkiId'")
        #self.connection = psycopg2.connect("host='localhost' port=5432 dbname='potential_fitting' user='USER' password='password'")
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

    def create(self):
        """
        Creates all the tables, procedures, and sequences in the database.

        Args:
            None.

        Returns:
            None.
        """

        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "init.sql")) as sql_initilizer:
            sql_script = sql_initilizer.read()

            self.cursor.execute(sql_script)

    def annihilate(self, confirm = "no way"):
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
            self.cursor.execute("TRUNCATE atom_info, fragment_contents, fragment_info, log_files, model_info, molecule_contents, molecule_info, molecule_list, molecule_properties, optimized_geometries, tags")
        else:
            print("annihilate failed. specify confirm = \"confirm\" if deletion of all content in the database is desired.")

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

        Args:
            command         - The command to run.
            params          - Parameters for the command.

        Returns:
            None.
        """

        self.cursor.execute(
            "DO $$" + 
            "   BEGIN " + 
            command + 
            "END $$",
            params)

    #TODO: exists, insert, select, and update will be removed in the future.

    def exists(self, table, **property_value_pairs):
        if len(property_value_pairs) < 1:
            # must specify at least one property
            raise InvalidValueError("property_value_pairs", property_value_pairs, "Must specify at least one pair.")
        properties = tuple(property_value_pairs)
        properties_string = properties[0] + "=%s"
        for property in properties[1:]:
            properties_string += " AND "
            properties_string += property + "=%s"

        values = tuple(property_value_pairs.values())
        self.commands.append("")
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
            raise InvalidValueError("property_value_pairs", property_value_pairs, "Must specify at least one pair.")
        if len(values) < 1:
            # must specify at least one value
            raise InvalidValueError("values", values, "Must specify at least one pair.")

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

    def update(self, table, values, set_property_value_pairs, where_property_value_pairs):
        if len(set_property_value_pairs) < 1:
            # must specify at least one property
            raise InvalidValueError("set_property_value_pairs", set_property_value_pairs, "Must specify at least one pair.")
        if len(where_property_value_pairs) < 1:
            # must specify at least one property
            raise InvalidValueError("where_property_value_pairs", where_property_value_pairs, "Must specify at least one pair.")

        raise NotImplementedError
        
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

    def add_model_info(self, method, basis, cp):\

        command_string = "PERFORM add_model_info(%s, %s, %s);"

        params = [method, basis, cp]

        return command_string, params

    def add_molecule(self, molecule):
        # should always be at COM and Principle axes before calling!!!!

        coordinates = []
        for fragment in molecule.get_fragments():
            for atom in fragment.get_atoms():
                coordinates.append(atom.get_x())
                coordinates.append(atom.get_y())
                coordinates.append(atom.get_z())


        command_string = "PERFORM add_molecule(%s, %s, ARRAY["
        params = (molecule.get_SHA1(), molecule.get_name())

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

            command_string += "construct_fragment(%s, %s, %s, %s, %s, %s)"
            if not frag_name == frag_names[-1]:
                command_string += ", "
            params += (frag_name, fragment.get_charge(), fragment.get_spin_multiplicity(), self.create_postgres_array(*symbols),
                    self.create_postgres_array(*symmetries), self.create_postgres_array(*counts))

        command_string += "], %s, %s);"
        params += (fragment_counts, coordinates)

        return command_string, params

    def add_calculations(self, molecule_list, method, basis, cp, *tags, optimized = False):

        command_string = ""
        params = []

        batch_count = 0

        for molecule in molecule_list:

            model_name ="{}/{}/{}".format(method, basis, cp)


            command_string += \
                "IF NOT EXISTS(SELECT mol_hash FROM molecule_properties WHERE mol_hash=%s AND model_name=%s) THEN "

            params += [molecule.get_SHA1(), model_name]

            new_command_string, new_params = self.add_molecule(molecule)
            command_string += new_command_string
            params += new_params

            new_command_string, new_params = self.add_model_info(method, basis, cp)
            command_string += new_command_string
            params += new_params

            for n in range(1, molecule.get_num_fragments() + 1):
                for fragment_indices in itertools.combinations(range(molecule.get_num_fragments()), n):
                    if cp and len(fragment_indices) != molecule.get_num_fragments():
                        command_string += "INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp) VALUES (%s, %s, %s, %s, %s, %s, %s, %s);"
                        params += [molecule.get_SHA1(), model_name, list(fragment_indices), [], [], "pending", [], True]

                    command_string += "INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp) VALUES (%s, %s, %s, %s, %s, %s, %s, %s);"
                    params += [molecule.get_SHA1(), model_name, list(fragment_indices), [], [], "pending", [], False]

            command_string += "INSERT INTO tags VALUES (%s, %s, %s);"
            params += [molecule.get_SHA1(), model_name, []]

            command_string += "END IF;"

            new_command_string, new_params = self.update_tags(molecule.get_SHA1(), model_name, *tags)
            command_string += new_command_string
            params += new_params
            
            if optimized:
                command_string += \
                    "IF NOT EXISTS(SELECT mol_hash FROM optimized_geometries WHERE mol_name=%s AND mol_hash=%s AND model_name=%s) THEN " + \
                    "   INSERT INTO optimized_geometries VALUES (%s, %s, %s);" + \
                    "END IF;"

                params += [molecule.get_name(), molecule.get_SHA1(), model_name, molecule.get_name(), molecule.get_SHA1(), model_name]

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
                self.execute(command_string, params)

    def update_tags(self, mol_hash, model_name, *tags):

        command_string = ""
        params = []

        for tag in tags:
            command_string += \
                "IF NOT EXISTS(SELECT mol_hash FROM tags WHERE mol_hash=%s AND model_name=%s AND %s=ANY(tag_names)) THEN " + \
                "   UPDATE tags SET tag_names=(%s || tag_names) WHERE mol_hash=%s AND model_name=%s;" + \
                "END IF;"

            params += [mol_hash, model_name, tag, self.create_postgres_array(tag), mol_hash, model_name]

        return command_string, params

    def build_empty_molecule(self, mol_name):
        """
        Returns a copy of the mol_name molecule from inside the database with all atom
        coordinates set to 0.

        Params:
            mol_name        - The name of the molecule to construct.

        Returns:
            a Molecule object with all coordinates set to 0.
        """

        molecule = Molecule()

        self.cursor.execute("SELECT frag_name, count, charge, spin FROM molecule_contents INNER JOIN fragment_info ON molecule_contents.frag_name=fragment_info.name where mol_name=%s", (mol_name,))

        molecule_contents = self.cursor.fetchall()

        for frag_name, frag_count, charge, spin in molecule_contents:
            for i in range(frag_count):
                fragment = Fragment(frag_name, charge, spin)

                fragment_contents = self.select("fragment_contents", True, "atom_symbol", "symmetry", "count", frag_name = frag_name)

                for atom_symbol, atom_symmetry, atom_count in fragment_contents:
                    for k in range(atom_count):
                        atom = Atom(atom_symbol, atom_symmetry, 0, 0, 0)

                        fragment.add_atom(atom)

                molecule.add_fragment(fragment)

        return molecule

    def get_all_calculations(self, client_name, *tags, calculations_to_do = sys.maxsize):

        # TODO: make tags work with this!

        while True:

            self.cursor.execute("SELECT molecule_list.mol_name from molecule_properties INNER JOIN molecule_list ON molecule_properties.mol_hash = molecule_list.mol_hash WHERE molecule_properties.status = %s", ("pending",))

            try:
                molecule_name = self.cursor.fetchone()[0]
            except TypeError:
                break

            self.cursor.execute("SELECT * FROM get_pending_calculations(%s, %s, %s, %s)", (molecule_name, client_name, self.create_postgres_array(*tags), min(self.batch_size, calculations_to_do)))

            calculations_to_do -= self.batch_size

            pending_calcs = self.cursor.fetchall()

            empty_molecule = self.build_empty_molecule(molecule_name)

            for atom_coordinates, model, frag_indices, use_cp in pending_calcs:
                molecule = copy.deepcopy(empty_molecule)

                for atom in molecule.get_atoms():
                    atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
                    atom_coordinates = atom_coordinates[3:]

                method = model[:model.index("/")]
                model = model[model.index("/") + 1:]
                basis = model[:model.index("/")]
                cp = model[model.index("/") + 1:]

                yield molecule, method, basis, cp, use_cp, frag_indices

            if calculations_to_do < 1:
                return

    def set_properties(self, calculation_results):
        # possibly check to make sure molecule is not morphed

        command_string = ""
        params = []

        batch_count = 0

        for molecule, method, basis, cp, use_cp, frag_indices, result, energy, log_text in calculation_results:
            model_name = method + "/" + basis + "/" + cp

            command_string += "PERFORM set_properties(%s, %s, %s, %s, %s, %s, %s);"
            params += [molecule.get_SHA1(), model_name, use_cp, self.create_postgres_array(*frag_indices), result, energy, log_text]

            batch_count += 1

            if batch_count == self.batch_size:
                self.execute(command_string, params)
                command_string = ""
                params = []
                batch_count = 0

        if batch_count != 0:
            self.execute(command_string, params)

    def get_1B_training_set(self, molecule_name, method, basis, cp, *tags):
        model_name ="{}/{}/{}".format(method, basis, cp)
        batch_offset = 0

        self.cursor.execute("SELECT * FROM count_entries(%s)", (molecule_name,))
        max_count = self.cursor.fetchone()[0]

        empty_molecule = self.build_empty_molecule(molecule_name)

        while True:
            self.cursor.execute("SELECT * FROM get_1B_training_set(%s, %s, %s, %s, %s)", (molecule_name, model_name, self.create_postgres_array(*tags), batch_offset, self.batch_size))
            training_set = self.cursor.fetchall()

            for atom_coordinates, energy in training_set:
                molecule = copy.deepcopy(empty_molecule)

                for atom in molecule.get_atoms():
                    atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
                    atom_coordinates = atom_coordinates[3:]

                yield molecule, energy

            batch_offset += self.batch_size

            if batch_offset > max_count:
                return;

    def get_2B_training_set(self, molecule_name, monomer1_name, monomer2_name, method, basis, cp, *tags):
        model_name ="{}/{}/{}".format(method, basis, cp)
        batch_offset = 0

        empty_molecule = self.build_empty_molecule(molecule_name)
        
        while True:
            self.cursor.execute("SELECT * FROM get_2B_training_set(%s, %s, %s, %s, %s, %s, %s)", (molecule_name, monomer1_name, monomer2_name, model_name, self.create_postgres_array(*tags), batch_offset, self.batch_size))
            training_set = self.cursor.fetchall()

            print("SET SIZE:", len(training_set))

            if training_set == []:
                return

            for atom_coordinates, binding_energy, interaction_energy, monomer1_energy, monomer2_energy in training_set:
                molecule = copy.deepcopy(empty_molecule)

                for atom in molecule.get_atoms():
                    atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
                    atom_coordinates = atom_coordinates[3:]

                yield molecule, binding_energy, interaction_energy, monomer1_energy, monomer2_energy

            batch_offset += self.batch_size

    def reset_all_calculations(self, *tags):
        self.cursor.execute("UPDATE molecule_properties SET status=%s WHERE status=%s", ("pending", "dispatched"))
        self.cursor.execute("UPDATE molecule_properties SET status=%s WHERE status=%s", ("pending", "complete"))

    def clean(self, *tags):
        self.cursor.execute("UPDATE molecule_properties SET status=%s WHERE status=%s", ("pending", "dispatched"))