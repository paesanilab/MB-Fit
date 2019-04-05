# external package imports
import itertools, psycopg2, numpy as np, copy, sys

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
            batch_size      - number of operations to perfrom on the database per roundtrip to the server.
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
        Simply returns self.

        Args:
            None.

        Returns:
            self.
        """

        return self

    def __exit__(self, type, value, traceback):
        """
        Saves and closes the database.

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
        Sets the number of operations to do on the database per server roundtrip.
        Larger numbers will be more efficient, but you should not exceed a couple thousand.

        Args:
            batch_size      - The number of operations to do per server roundtrip.

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
        Creates all the tables in the database, if they do not exist already.

        Also creates procedures and sequences.

        Args:
            None.

        Returns:
            None.
        """

        self.cursor.execute("""

            create schema public;

            comment on schema public is 'standard public schema';

            alter schema public owner to postgres;

            create type dup_result as
            (
                f1 integer,
                f2 text
            );

            alter type dup_result owner to ebullvul;

            create type empty_mol as
            (
                mol_name varchar,
                frag_names character varying[],
                frag_counts integer[],
                frag_charges integer[],
                frag_spins integer[]
            );

            alter type empty_mol owner to ebullvul;

            create table molecule_info
            (
                name varchar not null
                    constraint "molecule_Name_index"
                        primary key
            );

            alter table molecule_info owner to "paesanilab-pgadmins";

            create index molecule_info_name_uindex
                on molecule_info (name);

            create table molecule_list
            (
                mol_hash varchar not null
                    constraint molecule_list_hash_key
                        primary key,
                mol_name varchar not null
                    constraint molecule_list_molecule_name_fk
                        references molecule_info,
                atom_coordinates double precision[] not null
            );

            alter table molecule_list owner to "paesanilab-pgadmins";

            create unique index molecule_list_mol_hash_uindex
                on molecule_list (mol_hash);

            create index molecule_list_mol_name_index
                on molecule_list (mol_name);

            create table atom_info
            (
                atomic_symbol varchar not null
                    constraint atom_info_atomic_symbol_key
                        primary key
            );

            comment on table atom_info is 'Holds immutable information that is universal to all instances of a specific Atom.';

            alter table atom_info owner to "paesanilab-pgadmins";

            create unique index atom_info_atomic_symbol_uindex
                on atom_info (atomic_symbol);

            create table fragment_info
            (
                name varchar not null
                    constraint fragment_info_name_key
                        primary key,
                charge integer not null,
                spin integer not null
            );

            comment on table fragment_info is 'This table holds all immutable information shared by all fragments of a type';

            alter table fragment_info owner to "paesanilab-pgadmins";

            create unique index fragment_info_name_uindex
                on fragment_info (name);

            create table molecule_contents
            (
                mol_name varchar not null
                    constraint molecule_contents_molecule_info_name_fk
                        references molecule_info,
                frag_name varchar not null
                    constraint molecule_contents_fragment_info_name_fk
                        references fragment_info,
                count integer not null
            );

            alter table molecule_contents owner to "paesanilab-pgadmins";

            create index molecule_contents_mol_name_index
                on molecule_contents (mol_name);

            create table fragment_contents
            (
                frag_name varchar not null
                    constraint fragment_contents_fragment_info_name_fk
                        references fragment_info,
                atom_symbol varchar not null
                    constraint fragment_contents_atom_info_atomic_symbol_fk
                        references atom_info,
                count integer not null,
                symmetry varchar not null
            );

            comment on table fragment_contents is 'This table describes fragments found in fragment_info table as built by atoms in atom_info table.';

            alter table fragment_contents owner to "paesanilab-pgadmins";

            create index fragment_contents_frag_name_index
                on fragment_contents (frag_name);

            create table model_info
            (
                name varchar not null
                    constraint models_info_name_key
                        primary key
            );

            alter table model_info owner to "paesanilab-pgadmins";

            create unique index models_info_name_uindex
                on model_info (name);

            create table log_files
            (
                log_id serial not null
                    constraint log_files_log_id_key
                        primary key,
                start_time timestamp with time zone not null,
                end_time timestamp with time zone,
                log_text varchar,
                client_name varchar not null
            );

            alter table log_files owner to "paesanilab-pgadmins";

            create table molecule_properties
            (
                mol_hash varchar not null
                    constraint energies_list_molecule_list_mol_hash_fk
                        references molecule_list,
                model_name varchar not null
                    constraint energies_list_models_info_name_fk
                        references model_info,
                frag_indices integer[] not null,
                energies double precision[] not null,
                atomic_charges double precision[] not null,
                status varchar not null,
                past_log_ids integer[] not null,
                most_recent_log_id integer
                    constraint energies_list_log_files_log_id_fk
                        references log_files
            );

            alter table molecule_properties owner to "paesanilab-pgadmins";

            create index energies_list_mol_hash_model_name_index
                on molecule_properties (mol_hash, model_name);

            create unique index energies_list_most_recent_log_id_uindex
                on molecule_properties (most_recent_log_id);

            create index molecule_properties_status_index
                on molecule_properties (status);

            create index molecule_properties_mol_hash_index
                on molecule_properties (mol_hash);

            create index molecule_properties_frag_indices_index
                on molecule_properties (frag_indices);

            create index molecule_properties_mol_hash_model_name_frag_indices_index
                on molecule_properties (mol_hash, model_name, frag_indices);

            create unique index log_files_log_id_uindex
                on log_files (log_id);

            create table tags
            (
                mol_hash varchar not null
                    constraint tags_molecule_list_mol_hash_fk
                        references molecule_list,
                model_name varchar not null
                    constraint tags_models_info_name_fk
                        references model_info,
                tag_names character varying[] not null
            );

            alter table tags owner to "paesanilab-pgadmins";

            create unique index tags_mol_hash_model_name_uindex
                on tags (mol_hash, model_name);

            create table optimized_geometries
            (
                mol_name varchar not null
                    constraint optimized_geometries_molecule_info_name_fk
                        references molecule_info,
                mol_hash varchar not null
                    constraint optimized_geometries_molecule_list_mol_hash_fk
                        references molecule_list,
                model_name varchar
                    constraint optimized_geometries_models_info_name_fk
                        references model_info
            );

            alter table optimized_geometries owner to "paesanilab-pgadmins";

            create index optimized_geometries_mol_name_model_name_index
                on optimized_geometries (mol_name, model_name);

            create function get_empty_molecule(molecule_name character varying) returns empty_mol
                language plpgsql
            as $$
            DECLARE
                mol empty_mol;
                fragment_desc RECORD;
              BEGIN
              mol.mol_name := molecule_name;
              FOR fragment_desc IN SELECT frag_name, count, charge, spin
                  FROM molecule_contents INNER JOIN fragment_info ON molecule_contents.frag_name=fragment_info.name
                  where molecule_contents.mol_name=molecule_name LOOP

                  mol.frag_names := mol.frag_names || fragment_desc.frag_name;
                  mol.frag_counts := mol.frag_counts || fragment_desc.count;
                  mol.frag_charges := mol.frag_charges || fragment_desc.charge;
                  mol.frag_spins := mol.frag_spins || fragment_desc.spin;

                  --SELECT atom_symbol, symmetry, count FROM fragment_contents WHERE frag_name = fragment_desc.name INTO fragment_cont;

              END LOOP;
              RETURN mol;
              END;
            $$;

            alter function get_empty_molecule(varchar) owner to ebullvul;

            create function get_molecule(molecule_hash character varying) returns molecule_list
                language plpgsql
            as $$
            DECLARE
                mol empty_mol;
                fragment_desc RECORD;
              BEGIN
              mol.mol_name := molecule_name;
              FOR fragment_desc IN SELECT frag_name, count, charge, spin
                  FROM molecule_contents INNER JOIN fragment_info ON molecule_contents.frag_name=fragment_info.name
                  where molecule_contents.mol_name=molecule_name LOOP

                  mol.frag_names := mol.frag_names || fragment_desc.frag_name;
                  mol.frag_counts := mol.frag_counts || fragment_desc.count;
                  mol.frag_charges := mol.frag_charges || fragment_desc.charge;
                  mol.frag_spins := mol.frag_spins || fragment_desc.spin;

                  --SELECT atom_symbol, symmetry, count FROM fragment_contents WHERE frag_name = fragment_desc.name INTO fragment_cont;

              END LOOP;
              RETURN mol;
              END;
            $$;

            alter function get_molecule(varchar) owner to ebullvul;

            create function get_1b_training_set(molecule_name character varying, model character varying, input_tags character varying[], batch_offset integer, batch_size integer) returns TABLE(coords double precision[], energy double precision)
                language plpgsql
            as $$
            DECLARE
                    optimized_energy FLOAT = NULL;
                    hash VARCHAR;
                    mol_tags VARCHAR[];
                    mol_tag VARCHAR;
                    valid_tags BOOL;
                    optimized_properties molecule_properties;
                    mol_properties molecule_properties;
                  BEGIN

                    FOR hash IN SELECT mol_hash FROM optimized_geometries WHERE mol_name = molecule_name AND model_name = model
                    LOOP

                      valid_tags := FALSE;


                      SELECT tag_names FROM tags WHERE mol_hash = hash AND model_name = model LIMIT 1
                        INTO mol_tags;

                      FOREACH mol_tag IN ARRAY mol_tags
                      LOOP


                        IF (SELECT mol_tag = ANY(input_tags)) THEN
                          valid_tags := TRUE;
                        END IF;

                      END LOOP;


                      IF (SELECT valid_tags) THEN
                        IF (SELECT optimized_energy ISNULL) THEN
                          SELECT * FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0}'
                            INTO optimized_properties;
                          IF optimized_properties.status = 'complete' THEN
                            optimized_energy = optimized_properties.energies[1];
                          ELSE
                            RAISE EXCEPTION 'Optimized energy uncalculated in database.';
                          END IF;
                        ELSE
                          RAISE EXCEPTION 'Multiple optimized geometries in database.';
                        END IF;
                      END IF;

                    END LOOP;


                    IF (SELECT optimized_energy ISNULL) THEN
                      RAISE EXCEPTION 'No optimized energy in database.';
                    END IF;

                    FOR hash IN SELECT mol_hash FROM molecule_list WHERE mol_name = molecule_name OFFSET batch_offset LIMIT batch_size
                    LOOP
                      valid_tags := FALSE;


                      SELECT tag_names FROM tags WHERE mol_hash = hash AND model_name = model LIMIT 1
                        INTO mol_tags;

                      FOREACH mol_tag IN ARRAY mol_tags
                      LOOP


                        IF (SELECT mol_tag = ANY(input_tags)) THEN
                          valid_tags := TRUE;
                        END IF;

                      END LOOP;

                      IF (SELECT valid_tags) THEN
                        SELECT * FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0}'
                            INTO mol_properties;
                          IF mol_properties.status = 'complete' THEN
                            SELECT atom_coordinates  FROM molecule_list WHERE mol_hash = hash
                              INTO coords;
                            energy := mol_properties.energies[1] - optimized_energy;
                            RETURN NEXT;
                          END IF;
                      END IF;

                    END LOOP;
                  END;

            $$;

            alter function get_1b_training_set(varchar, varchar, character varying[], integer, integer) owner to ebullvul;

            create function set_properties(hash character varying, model character varying, indices integer[], energy double precision, log_txt character varying) returns void
                language plpgsql
            as $$
            DECLARE
                id INTEGER;

              BEGIN

                UPDATE molecule_properties SET energies = energies || energy, status='complete' WHERE mol_hash=hash AND model_name=model AND frag_indices=indices RETURNING most_recent_log_id
                  INTO id;
                UPDATE log_files SET end_time=clock_timestamp(), log_text=log_txt WHERE log_id=id;

              END;

            $$;

            alter function set_properties(varchar, varchar, integer[], double precision, varchar) owner to ebullvul;

            create function get_pending_calculations(molecule_name character varying, input_client_name character varying, input_tags character varying[], batch_size integer) returns TABLE(coords double precision[], model character varying, indices integer[])
                language plpgsql
            as $$
            DECLARE
                molecule_hash VARCHAR;
                id INTEGER;

              BEGIN

                FOR molecule_hash, model, indices IN SELECT molecule_properties.mol_hash, molecule_properties.model_name, molecule_properties.frag_indices FROM
                  molecule_properties INNER JOIN molecule_list ON molecule_properties.mol_hash = molecule_list.mol_hash WHERE
                  molecule_properties.status = 'pending' AND molecule_list.mol_name = molecule_name LIMIT batch_size
                LOOP

                  INSERT INTO log_files(start_time, client_name) VALUES (clock_timestamp(), input_client_name) RETURNING log_id
                    INTO id;

                  UPDATE molecule_properties SET status='dispatched', most_recent_log_id=id WHERE mol_hash = molecule_hash AND model_name = model AND frag_indices = indices;

                  SELECT atom_coordinates, mol_name from molecule_list WHERE mol_hash = molecule_hash
                    INTO coords;

                  RETURN NEXT;

                END LOOP;

              END;

            $$;

            alter function get_pending_calculations(varchar, varchar, character varying[], integer) owner to ebullvul;

            create function get_2b_training_set(molecule_name character varying, monomer1_name character varying, monomer2_name character varying, model character varying, input_tags character varying[], batch_offset integer, batch_size integer) returns TABLE(coords double precision[], binding_energy double precision, interaction_energy double precision, monomer1_deformation_energy double precision, monomer2_deformation_energy double precision)
                language plpgsql
            as $$
            DECLARE
                    optimized_monomer1_energy FLOAT = NULL;
                    optimized_monomer2_energy FLOAT = NULL;
                    hash VARCHAR;
                    mol_tags VARCHAR[];
                    mol_tag VARCHAR;
                    valid_tags BOOL;
                    optimized_properties molecule_properties;
                    dimer_status VARCHAR;
                    dimer_energies FLOAT[];
                    monomer1_energies FLOAT[];
                    monomer2_energies FLOAT[];
                    coordinates FLOAT[];
                  BEGIN

                    FOR hash IN SELECT mol_hash FROM optimized_geometries WHERE mol_name = monomer1_name AND model_name = model
                    LOOP

                      valid_tags := FALSE;


                      SELECT tag_names FROM tags WHERE mol_hash = hash AND model_name = model LIMIT 1
                        INTO mol_tags;

                      FOREACH mol_tag IN ARRAY mol_tags
                      LOOP


                        IF (SELECT mol_tag = ANY(input_tags)) THEN
                          valid_tags := TRUE;
                        END IF;

                      END LOOP;


                      IF (SELECT valid_tags) THEN
                        IF (SELECT optimized_monomer1_energy ISNULL) THEN
                          SELECT * FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0}'
                            INTO optimized_properties;
                          IF optimized_properties.status = 'complete' THEN
                            optimized_monomer1_energy = optimized_properties.energies[1];
                          ELSE
                            RAISE EXCEPTION 'Optimized energy uncalculated in database.';
                          END IF;
                        ELSE
                          RAISE EXCEPTION 'Multiple optimized geometries in database.';
                        END IF;
                      END IF;

                    END LOOP;

                    IF (SELECT optimized_monomer1_energy ISNULL) THEN
                      RAISE EXCEPTION 'No monomer1 optimized energy in database.';
                    END IF;

                    FOR hash IN SELECT mol_hash FROM optimized_geometries WHERE mol_name = monomer2_name AND model_name = model
                    LOOP

                      valid_tags := FALSE;


                      SELECT tag_names FROM tags WHERE mol_hash = hash AND model_name = model LIMIT 1
                        INTO mol_tags;

                      FOREACH mol_tag IN ARRAY mol_tags
                      LOOP


                        IF (SELECT mol_tag = ANY(input_tags)) THEN
                          valid_tags := TRUE;
                        END IF;

                      END LOOP;


                      IF (SELECT valid_tags) THEN
                        IF (SELECT optimized_monomer2_energy ISNULL) THEN
                          SELECT * FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0}'
                            INTO optimized_properties;
                          IF optimized_properties.status = 'complete' THEN
                            optimized_monomer2_energy = optimized_properties.energies[1];
                          ELSE
                            RAISE EXCEPTION 'Optimized energy uncalculated in database.';
                          END IF;
                        ELSE
                          RAISE EXCEPTION 'Multiple optimized geometries in database.';
                        END IF;
                      END IF;

                    END LOOP;

                    IF (SELECT optimized_monomer2_energy ISNULL) THEN
                      RAISE EXCEPTION 'No monomer2 optimized energy in database.';
                    END IF;

                    FOR hash IN SELECT mol_hash FROM molecule_list WHERE mol_name = molecule_name OFFSET batch_offset LIMIT batch_size
                    LOOP
                      valid_tags := FALSE;

                      SELECT tag_names FROM tags WHERE mol_hash = hash AND model_name = model LIMIT 1
                        INTO mol_tags;

                      FOREACH mol_tag IN ARRAY mol_tags
                      LOOP


                        IF (SELECT mol_tag = ANY(input_tags)) THEN
                          valid_tags := TRUE;
                        END IF;

                      END LOOP;

                      IF (SELECT valid_tags) THEN
                        SELECT energies, status FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0, 1}'
                            INTO dimer_energies, dimer_status;

                        IF dimer_status = 'complete' THEN
                          SELECT atom_coordinates  FROM molecule_list WHERE mol_hash = hash
                            INTO coords;

                          SELECT energies FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0}'
                            INTO monomer1_energies;
                          SELECT energies FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{1}'
                            INTO monomer2_energies;

                          interaction_energy := dimer_energies[1] - monomer1_energies[1] - monomer2_energies[1];
                          monomer1_deformation_energy := monomer1_energies[1] - optimized_monomer1_energy;
                          monomer2_deformation_energy := monomer2_energies[1] - optimized_monomer2_energy;
                          binding_energy := interaction_energy - monomer1_deformation_energy - monomer2_deformation_energy;
                          RETURN NEXT;
                        END IF;
                      END IF;

                    END LOOP;
                  END;

            $$;

            alter function get_2b_training_set(varchar, varchar, varchar, varchar, character varying[], integer, integer) owner to ebullvul;

        """)

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
        One param per '%s'

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
        Creates a postgres array with an arbitrary number of values

        Args:
            *values: Set of abritrary number of values which firm the lements of the postgres array.

        Returns:
            The Postgres Array as a python string 
        """

        return "{" + ",".join([str(i) for i in values]) + "}"

    def add_atom_info(self, atom):
        """
        Adds a single atom type's info to the atom_info table if it does not already exist.

        Args:
            atom: The atom to insert (Object of the Atom class)

        Returns:
            command_string: The postgres command which inserts the atom into the atom_info table as a Python String.
            params: A 2 element list containing the name of the atom twice.   

        """
        command_string = \
            "IF NOT EXISTS(SELECT atomic_symbol FROM atom_info WHERE atomic_symbol=%s) THEN " + \
            "    INSERT INTO atom_info VALUES (%s);" + \
            "END IF;"

        params = [atom.get_name(), atom.get_name()]

        return command_string, params
            
    def add_fragment_info(self, fragment):

        # use transactions

        """
        Adds a single fragment type's info to the fragment_info table if it does not already exist.

        Args:
            fragment: The fragment to insert (Object of the Fragment class)

        Returns:
            command_string: The postgres command which inserts the fragment into the fragment_contents table as a Python String.
            params: A list of the name of the fragment, the atom symbol of tthe first symmetry pair, the count associated with the fragment, and the symmetry symbol.  

        """

        command_string = ""
        params = []

        command_string += \
            "IF NOT EXISTS(SELECT name FROM fragment_info WHERE name=%s AND charge=%s AND spin=%s) THEN " + \
            "    INSERT INTO fragment_info VALUES(%s, %s, %s);"

        params += [fragment.get_name(), fragment.get_charge(), fragment.get_spin_multiplicity(), fragment.get_name(), fragment.get_charge(), fragment.get_spin_multiplicity()]

        for atom in fragment.get_atoms():
            new_command_string, new_params = self.add_atom_info(atom)
            command_string += new_command_string
            params += new_params

        # updates fragment_contents
        atoms = [[atom.get_name(), atom.get_symmetry_class()] for atom in fragment.get_atoms()]

        symbol_symmetry_pairs, counts = np.unique(atoms, return_counts = True, axis = 0)
        counts = [int(i) for i in counts]

        for symbol_symmetry_pair, count in zip(symbol_symmetry_pairs, counts):
            atomic_symbol = symbol_symmetry_pair[0]
            symmetry_symbol = symbol_symmetry_pair[1]
            command_string += "    INSERT INTO fragment_contents VALUES (%s, %s, %s, %s);"
            params += [fragment.get_name(), atomic_symbol, count, symmetry_symbol]

        command_string += "END IF;"

        return command_string, params

    def add_molecule_info(self, molecule):

        """
        Adds a single molecule type's info to the molecule_info table if it does not already exist
        """

        command_string = ""
        params = []

        command_string += \
            "IF NOT EXISTS(SELECT name FROM molecule_info WHERE name=%s) THEN " + \
            "   INSERT INTO molecule_info VALUES(%s);"
        params += [molecule.get_name(), molecule.get_name()]

        for fragment in molecule.get_fragments():
            new_command_string, new_params = self.add_fragment_info(fragment)
            command_string += new_command_string
            params += new_params



        fragments, counts = np.unique([fragment.get_name() for fragment in molecule.get_fragments()], return_counts = True)
        counts = [int(i) for i in counts]

        for fragment, count in zip(fragments, counts):
            command_string += "    INSERT INTO molecule_contents VALUES (%s, %s, %s);"
            params += [molecule.get_name(), fragment, count]

        command_string += "END IF;"

        return command_string, params

    def add_model_info(self, method, basis, cp):
        model_name ="{}/{}/{}".format(method, basis, cp)

        command_string = \
            "IF NOT EXISTS(SELECT name FROM model_info WHERE name=%s) THEN" + \
            "   INSERT INTO model_info VALUES (%s);" + \
            "END IF;"

        params = [model_name, model_name]

        return command_string, params

    def add_molecule(self, molecule):
        # should always be at COM and Principle axes before calling!!!!

        command_string = ""
        params = []

        command_string += \
            "IF NOT EXISTS(SELECT mol_hash FROM molecule_list WHERE mol_hash=%s) THEN "

        params += [molecule.get_SHA1()]

        new_command_string, new_params = self.add_molecule_info(molecule)
        command_string += new_command_string
        params += new_params



        coordinates = []
        for fragment in molecule.get_fragments():
            for atom in fragment.get_atoms():
                coordinates.append(atom.get_x())
                coordinates.append(atom.get_y())
                coordinates.append(atom.get_z())

        command_string += \
            "   INSERT INTO molecule_list VALUES(%s, %s, %s);" + \
            "END IF;"

        params += [molecule.get_SHA1(), molecule.get_name(), coordinates]

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
                    command_string += "INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids) VALUES (%s, %s, %s, %s, %s, %s, %s);"
                    params += [molecule.get_SHA1(), model_name, list(fragment_indices), [], [], "pending", []]

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

    def get_molecule(self, mol_hash, molecule = None):

        mol_name, atom_coordinates = self.select("molecule_list", False, "mol_name", "atom_coordinates", mol_hash = mol_hash)

        if molecule is None:
            molecule = self.build_empty_molecule(mol_name)

        for atom in molecule.get_atoms():
            atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
            atom_coordinates = atom_coordinates[3:]

        return molecule

    def build_empty_molecule(self, mol_name):
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

            for atom_coordinates, model, frag_indices in pending_calcs:
                molecule = copy.deepcopy(empty_molecule)

                for atom in molecule.get_atoms():
                    atom.set_xyz(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2])
                    atom_coordinates = atom_coordinates[3:]

                method = model[:model.index("/")]
                model = model[model.index("/") + 1:]
                basis = model[:model.index("/")]
                cp = model[model.index("/") + 1:]

                yield molecule, method, basis, cp, frag_indices

            if calculations_to_do < 1:
                return

    def set_properties(self, calculation_results):
        # possibly check to make sure molecule is not morphed

        command_string = ""
        params = []

        batch_count = 0

        for molecule, method, basis, cp, frag_indices, result, energy, log_text in calculation_results:
            model_name = method + "/" + basis + "/" + cp

            command_string += "PERFORM set_properties(%s, %s, %s, %s, %s, %s);"
            params += [molecule.get_SHA1(), model_name, self.create_postgres_array(*frag_indices), result, energy, log_text]

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