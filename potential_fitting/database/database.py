# external package imports
import itertools, psycopg2, numpy as np, copy, os

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

    def get_all_calculations(self, client_name, *tags):

        # TODO: make tags work with this!

        while True:

            self.cursor.execute("SELECT molecule_list.mol_name from molecule_properties INNER JOIN molecule_list ON molecule_properties.mol_hash = molecule_list.mol_hash WHERE molecule_properties.status = %s", ("pending",))

            try:
                molecule_name = self.cursor.fetchone()[0]
            except TypeError:
                break

            self.cursor.execute("SELECT * FROM get_pending_calculations(%s, %s, %s, %s)", (molecule_name, client_name, self.create_postgres_array(*tags), self.batch_size))

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
