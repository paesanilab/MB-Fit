-- we don't know how to generate root <with-no-name> (class Root) :(
create type empty_mol as
(
	mol_name varchar,
	frag_names character varying[],
	frag_counts integer[],
	frag_charges integer[],
	frag_spins integer[],
	frag_smiles character varying[]
);

create type fragment as
(
	name varchar,
	charge integer,
	spin integer,
	atomic_symbols character varying[],
	symmetries character varying[],
	counts integer[],
	smile varchar
);

create type status_enum as enum ('pending', 'dispatched', 'complete', 'failed');

create table molecule_info
(
	name varchar not null
		constraint "molecule_Name_index"
			primary key
);

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

create unique index atom_info_atomic_symbol_uindex
	on atom_info (atomic_symbol);

create table fragment_info
(
	name varchar not null
		constraint fragment_info_name_key
			primary key,
	charge integer not null,
	spin integer not null,
	smile varchar default 'no smile specified'::character varying not null
);

comment on table fragment_info is 'This table holds all immutable information shared by all fragments of a type';

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

create index fragment_contents_frag_name_index
	on fragment_contents (frag_name);

create table model_info
(
	name varchar not null
		constraint models_info_name_key
			primary key
);

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

create unique index log_files_log_id_uindex
	on log_files (log_id);

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
	status status_enum not null,
	past_log_ids integer[] not null,
	most_recent_log_id integer
		constraint energies_list_log_files_log_id_fk
			references log_files,
	use_cp boolean not null
);

create index energies_list_mol_hash_model_name_index
	on molecule_properties (mol_hash, model_name);

create unique index energies_list_most_recent_log_id_uindex
	on molecule_properties (most_recent_log_id);

create index molecule_properties_mol_hash_index
	on molecule_properties (mol_hash);

create index molecule_properties_frag_indices_index
	on molecule_properties (frag_indices);

create index molecule_properties_mol_hash_model_name_frag_indices_index
	on molecule_properties (mol_hash, model_name, frag_indices);

create index molecule_properties_status_index
	on molecule_properties (status);

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

create index optimized_geometries_mol_name_model_name_index
	on optimized_geometries (mol_name, model_name);

create table pending_calculations
(
	mol_hash varchar not null
		constraint pending_calculations_molecule_list_mol_hash_fk
			references molecule_list,
	model_name varchar not null,
	frag_indices integer[] not null,
	use_cp boolean not null
);

create index pending_calculations_mol_hash_model_name_frag_indices_use_cp_in
	on pending_calculations (mol_hash, model_name, frag_indices, use_cp);

create index pending_calculations_mol_hash_model_name_index
	on pending_calculations (mol_hash, model_name);

create table training_sets
(
	admins integer[] not null,
	tag_name varchar not null
		constraint training_sets_pk
			primary key,
	read_users integer[] not null,
	write_users integer[] not null
);

create unique index training_sets_tag_name_uindex
	on training_sets (tag_name);

create function get_molecule(molecule_hash character varying) returns molecule_list
	security definer
	SET search_path=public, pg_temp
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
      mol.frag_SMILES := mol.frag_SMILES || fragment_desc.SMILE;

      --SELECT atom_symbol, symmetry, count FROM fragment_contents WHERE frag_name = fragment_desc.name INTO fragment_cont;

  END LOOP;
  RETURN mol;
  END;
$$;

create function count_entries(molecule_name character varying) returns integer
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    count integer;
  BEGIN
    SELECT COUNT(*) FROM molecule_list WHERE mol_name=molecule_name
      into count;
    RETURN count;
  END;

$$;

create function add_atom_info(atomic_symbol character varying) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE

  BEGIN
    IF NOT EXISTS(SELECT atom_info.atomic_symbol FROM atom_info WHERE atom_info.atomic_symbol=add_atom_info.atomic_symbol) THEN
      INSERT INTO atom_info VALUES (add_atom_info.atomic_symbol);
      RETURN True;
    ELSE
      RETURN False;
    END IF;
  END;

$$;

create function add_molecule_info(name character varying, fragments fragment[], counts integer[]) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  index integer;

  BEGIN
    IF NOT EXISTS(SELECT molecule_info.name FROM molecule_info WHERE molecule_info.name=add_molecule_info.name) THEN
      INSERT INTO molecule_info VALUES(add_molecule_info.name);

      FOR index IN 1..array_length(fragments, 1) LOOP
        PERFORM add_fragment_info(fragments[index].name, fragments[index].charge, fragments[index].spin, fragments[index].smile, fragments[index].atomic_symbols, fragments[index].symmetries, fragments[index].counts);
        INSERT INTO molecule_contents VALUES (add_molecule_info.name, fragments[index].name, counts[index]);
      END LOOP;

      RETURN True;
    ELSE
      RETURN False;
    END IF;
  END;

$$;

create function add_model_info(method character varying, basis character varying, cp boolean) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  model_name varchar;
  BEGIN
    model_name := concat(method, '/', basis, '/');
    IF cp THEN
      model_name := concat(model_name, 'True');
    ELSE
      model_name := concat(model_name, 'False');
    end if;

    IF NOT EXISTS(SELECT name FROM model_info WHERE name=model_name) THEN
      INSERT INTO model_info VALUES (model_name);
      RETURN True;
    ELSE
      RETURN False;
    END IF;
  END;

$$;

create function add_molecule(hash character varying, name character varying, fragments fragment[], counts integer[], coordinates double precision[]) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE

  BEGIN
    IF NOT EXISTS(SELECT mol_hash from molecule_list WHERE mol_hash=hash) THEN
      PERFORM add_molecule_info(name, fragments, counts);
      INSERT INTO molecule_list VALUES(hash, name, coordinates);
      RETURN True;
    ELSE
      RETURN False;
    END IF;
  END;

$$;

create function combinations(arr integer[]) returns TABLE(perm integer[], l integer)
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  len INTEGER;
  x INTEGER[];
  i INTEGER;
  BEGIN
    SELECT * FROM array_length(arr, 1) INTO len;
    IF len = 1 THEN
      perm := arr;
      l := 1;
      RETURN NEXT;
    ELSE
      perm := '{}';
      perm := (perm || arr[1]);
      l := 1;
      RETURN NEXT;

      FOR x IN SELECT * FROM combinations(arr[2:]) LOOP
        perm := x;
        l := array_length(perm, 1);
        RETURN NEXT;
      END LOOP;

      FOR x IN SELECT * FROM combinations(arr[2:]) LOOP
        perm := (arr[1] || x);
        l := array_length(perm, 1);
        RETURN NEXT;
      END LOOP;

    END IF;
  END;

$$;

create function reset_dispatched(ts character varying[]) returns integer
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    hash VARCHAR;
    model VARCHAR;
    frags INTEGER[];
    cp BOOLEAN;
    counter INTEGER;
  BEGIN
    counter := 0;
    FOR hash, model, frags, cp IN SELECT molecule_properties.mol_hash, molecule_properties.model_name, molecule_properties.frag_indices, molecule_properties.use_cp FROM molecule_properties INNER JOIN tags ON molecule_properties.mol_hash = tags.mol_hash AND molecule_properties.model_name = tags.model_name WHERE status = 'dispatched' and ts && tags.tag_names
    LOOP
      UPDATE molecule_properties SET status = 'pending' WHERE mol_hash = hash AND model_name = model AND frag_indices = frags;
      INSERT INTO pending_calculations VALUES (hash, model, frags, cp);
      counter = counter + 1;
    end loop;

    RETURN counter;
  END;

$$;

create function reset_failed(ts character varying[]) returns integer
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    hash VARCHAR;
    model VARCHAR;
    frags INTEGER[];
    cp BOOLEAN;
    counter INTEGER;
  BEGIN
    counter := 0;
    FOR hash, model, frags, cp IN SELECT molecule_properties.mol_hash, molecule_properties.model_name, molecule_properties.frag_indices, molecule_properties.use_cp FROM molecule_properties INNER JOIN tags ON molecule_properties.mol_hash = tags.mol_hash AND molecule_properties.model_name = tags.model_name WHERE status = 'failed' and ts && tags.tag_names
    LOOP
      UPDATE molecule_properties SET status = 'pending' WHERE mol_hash = hash AND model_name = model AND frag_indices = frags;
      INSERT INTO pending_calculations VALUES (hash, model, frags, cp);
      counter = counter + 1;
    end loop;

    RETURN counter;
  END;

$$;

create function add_fragment_info(name character varying, charge integer, spin integer, smile character varying, atomic_symbols character varying[], symmetries character varying[], counts integer[]) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  atomic_symbol character varying;
  index integer;
  BEGIN
    IF NOT EXISTS(SELECT fragment_info.name FROM fragment_info WHERE fragment_info.name=add_fragment_info.name AND fragment_info.charge=add_fragment_info.charge AND fragment_info.spin=add_fragment_info.spin AND fragment_info.smile=add_fragment_info.SMILE) THEN
      INSERT INTO fragment_info VALUES(add_fragment_info.name, add_fragment_info.charge, add_fragment_info.spin, add_fragment_info.SMILE);

      FOREACH atomic_symbol IN ARRAY atomic_symbols LOOP
        PERFORM add_atom_info(atomic_symbol);
      END LOOP;

      FOR index IN 1..array_length(atomic_symbols, 1) LOOP
        INSERT INTO fragment_contents VALUES (name, atomic_symbols[index], counts[index], symmetries[index]);
      END LOOP;

      RETURN True;
    ELSE
      RETURN False;
    END IF;
  END;

$$;

create function construct_fragment(name character varying, charge integer, spin integer, smile character varying, atomic_symbols character varying[], symmetries character varying[], counts integer[]) returns fragment
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  frag fragment;

  BEGIN
    frag.name = name;
    frag.charge = charge;
    frag.spin = spin;
    frag.atomic_symbols = atomic_symbols;
    frag.symmetries = symmetries;
    frag.counts = counts;
    frag.smile = SMILE;
    RETURN frag;
  END;

$$;

create function get_pending_calculations(molecule_name character varying, input_client_name character varying, input_tags character varying[], batch_size integer) returns TABLE(coords double precision[], model character varying, indices integer[], use_cp boolean)
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    molecule_hash VARCHAR;
    id INTEGER;

  BEGIN

    FOR molecule_hash, model, indices, get_pending_calculations.use_cp IN SELECT pending_calculations.mol_hash, pending_calculations.model_name, pending_calculations.frag_indices, pending_calculations.use_cp FROM
      pending_calculations INNER JOIN molecule_list ON pending_calculations.mol_hash = molecule_list.mol_hash
      INNER JOIN tags ON pending_calculations.mol_hash = tags.mol_hash AND pending_calculations.model_name = tags.model_name WHERE
      molecule_list.mol_name = molecule_name AND tags.tag_names && input_tags LIMIT batch_size
    LOOP

      INSERT INTO log_files(start_time, client_name) VALUES (clock_timestamp(), input_client_name) RETURNING log_id
        INTO id;

      UPDATE molecule_properties SET status='dispatched', most_recent_log_id=id WHERE mol_hash = molecule_hash AND model_name = model AND frag_indices = indices AND molecule_properties.use_cp = get_pending_calculations.use_cp;

      DELETE FROM pending_calculations WHERE mol_hash = molecule_hash AND model_name = model AND frag_indices = indices AND pending_calculations.use_cp = get_pending_calculations.use_cp;

      SELECT atom_coordinates, mol_name from molecule_list WHERE mol_hash = molecule_hash
        INTO coords;

      RETURN NEXT;

    END LOOP;

  END;

$$;

create function delete_atom_info(atomic_symbol character varying) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE

  BEGIN
    IF NOT EXISTS(SELECT atom_symbol FROM fragment_contents WHERE atom_symbol = atomic_symbol) THEN
      DELETE FROM atom_info WHERE atom_info.atomic_symbol = delete_atom_info.atomic_symbol;
      RETURN TRUE;
    end if;
    RETURN FALSE;
  END;

$$;

create function delete_fragment_info(name character varying) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  atomic_symbol character varying;
  atom_symmetry character varying;
  BEGIN
    IF NOT EXISTS(SELECT frag_name FROM molecule_contents WHERE frag_name = name) THEN
      for atomic_symbol, atom_symmetry in SELECT atom_symbol, symmetry FROM fragment_contents WHERE frag_name = name LOOP
        DELETE FROM fragment_contents WHERE frag_name = name AND atomic_symbol = atomic_symbol AND symmetry = atom_symmetry;
        PERFORM delete_atom_info(atomic_symbol);
      end loop;

      DELETE FROM fragment_info WHERE fragment_info.name = delete_fragment_info.name;

      RETURN TRUE;
    end if;
    RETURN FALSE;
  END;

$$;

create function delete_molecule_info(name character varying) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  fragment_name varchar;

  BEGIN

    IF NOT EXISTS(SELECT mol_name FROM molecule_list WHERE mol_name=name) THEN
      for fragment_name in SELECT frag_name FROM molecule_contents WHERE mol_name = name LOOP
        DELETE FROM molecule_contents WHERE mol_name = name AND frag_name = fragment_name;
        PERFORM delete_fragment_info(fragment_name);
      end loop;

      DELETE FROM molecule_info WHERE molecule_info.name = delete_molecule_info.name;

      RETURN TRUE;
    end if;

    RETURN FALSE;
  END;

$$;

create function delete_molecule(hash character varying, name character varying) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE

  BEGIN

    IF NOT EXISTS(SELECT mol_hash FROM molecule_properties WHERE mol_hash=hash) THEN

      DELETE FROM molecule_list WHERE mol_hash=hash;

      PERFORM delete_molecule_info(name);

      RETURN TRUE;
    end if;

    RETURN FALSE;
  END;

$$;

create function delete_model_info(method character varying, basis character varying, cp boolean) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  model varchar;
  BEGIN
    model := concat(method, '/', basis, '/');
    IF cp THEN
      model := concat(model, 'True');
    ELSE
      model := concat(model, 'False');
    end if;

    IF NOT EXISTS(SELECT model_name FROM molecule_properties WHERE model_name = model) THEN
      DELETE FROM model_info WHERE name=model;

      RETURN TRUE;
    end if;

    RETURN FALSE;
  END;

$$;

create function delete_all_calculations(molecule_name character varying, method character varying, basis character varying, cp boolean, input_tags character varying[], delete_complete_calculations boolean) returns integer
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  model   varchar;
  hash    varchar;
  counter int;
BEGIN

  model := concat(method, '/', basis, '/');
  IF cp THEN
    model := concat(model, 'True');
  ELSE
    model := concat(model, 'False');
  end if;

  counter := 0;

  for hash in SELECT molecule_list.mol_hash
                     FROM tags
                            INNER JOIN molecule_list ON tags.mol_hash = molecule_list.mol_hash
                     WHERE molecule_list.mol_name = molecule_name
                       AND model_name = model
                       AND tag_names && input_tags
    LOOP

      IF delete_calculation(hash, molecule_name, method, basis, cp, input_tags, delete_complete_calculations) THEN
        counter = counter + 1;
      end if;

    end loop;

  return counter;
END;

$$;

create function training_set_exists(input_tag character varying) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  exists boolean;
      BEGIN

        SELECT EXISTS(SELECT tag_name FROM training_sets WHERE tag_name = input_tag)
          into exists;

        RETURN exists;
      END;

$$;

create function get_user_id(username character varying) returns integer
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  user_id int;
      BEGIN

        SELECT usesysid FROM pg_user
            WHERE usename=username
            INTO user_id;

        RETURN user_id;
      END;

$$;

create function has_read_privilege(input_tag character varying) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  privileged boolean;
      BEGIN

        SELECT get_user_id() = ANY(read_users) FROM training_sets WHERE tag_name = input_tag
          into privileged;

        RETURN privileged;
      END;

$$;

create function has_write_privilege(input_tag character varying) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  privileged boolean;
      BEGIN

        SELECT get_user_id() = ANY(write_users) FROM training_sets WHERE tag_name = input_tag
          into privileged;

        RETURN privileged;
      END;

$$;

create function has_admin_privilege(input_tag character varying) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  privileged boolean;
      BEGIN

        SELECT get_user_id() = ANY(admins) FROM training_sets WHERE tag_name = input_tag
          into privileged;

        RETURN privileged;
      END;

$$;

create function get_empty_molecule(molecule_name character varying) returns TABLE(f_name character varying, f_charge integer, f_spin integer, a_atomic_symbols character varying[], a_symmetries character varying[], a_counts integer[], f_smile character varying, f_count integer)
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    atom varchar;
    atom_symmetry varchar;
    num_atoms integer;
  BEGIN
    FOR f_name, f_count
        IN SELECT frag_name, count
        FROM molecule_contents
        where mol_name=molecule_name
    LOOP
      SELECT charge, spin, SMILE
          FROM fragment_info
          WHERE name = f_name
      INTO f_charge, f_spin, f_SMILE;

      a_atomic_symbols := '{}';
      a_symmetries := '{}';
      a_counts := '{}';

      FOR atom, atom_symmetry, num_atoms
          IN SELECT atom_symbol, symmetry, count
          FROM fragment_contents
          WHERE frag_name=f_name
      LOOP

        a_atomic_symbols := a_atomic_symbols || atom;
        a_symmetries := a_symmetries || atom_symmetry;
        a_counts := a_counts || num_atoms;

      END LOOP;

      RETURN NEXT;

    END LOOP;
  END;
$$;

create function annihilate() returns void
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE

  BEGIN
    TRUNCATE atom_info, fragment_contents, fragment_info, log_files, model_info, molecule_contents, molecule_info, molecule_list, molecule_properties, pending_calculations, optimized_geometries, tags, training_sets;
  END;

$$;

create function add_calculation(hash character varying, name character varying, fragments fragment[], counts integer[], coordinates double precision[], method character varying, basis character varying, cp boolean, tags character varying[], optimized boolean) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  model varchar;
  indices INTEGER[];
  all_frags INTEGER[];
  i INTEGER;
  x integer;
  z integer;
  new_config BOOLEAN;
  tag_name VARCHAR;
  BEGIN

    model := concat(method, '/', basis, '/');
    IF cp THEN
      model := concat(model, 'True');
    ELSE
      model := concat(model, 'False');
    end if;

    -- Check to Make sure the current user has write priveleges on all training sets for this geometry

    FOREACH tag_name IN ARRAY tags
    LOOP
      IF NOT training_set_exists(tag_name)
      THEN
        -- If the training set does not exist, create it giving current user all privileges
        INSERT INTO training_sets (admins, tag_name, read_users, write_users)
            VALUES (ARRAY[get_user_id()], tag_name, ARRAY[get_user_id()], ARRAY[get_user_id()]);
      ELSIF NOT has_write_privilege(tag_name)
      THEN
        -- If training set does exist and current user doesn't have write privileges, Error.
        raise EXCEPTION 'User %% does not have write privileges on training set %%', session_user, tag_name;
      END IF;
    END LOOP;

    IF NOT EXISTS(SELECT mol_hash FROM molecule_properties WHERE mol_hash=hash AND model_name=model) THEN
      PERFORM add_molecule(hash, name, fragments, counts, coordinates);
      PERFORM add_model_info(method, basis, cp);


      i = 0;
      FOREACH x IN ARRAY counts LOOP
        FOR z IN 1..x LOOP
          all_frags := all_frags || i;
          i := i + 1;
        END LOOP;
      END LOOP;

      FOR indices IN SELECT perm FROM combinations(all_frags) LOOP
        IF cp = True
        THEN

          IF array_length(indices, 1) = array_length(all_frags, 1)
          THEN
            INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
                VALUES (hash, model, indices, '{}', '{}', 'pending', '{}', False);

            INSERT INTO pending_calculations (mol_hash, model_name, frag_indices, use_cp)
                VALUES (hash, model, indices, False);
          ELSEIF array_length(indices, 1) = 1
          THEN
            INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
                VALUES (hash, model, indices, '{}', '{}', 'pending', '{}', True);

            INSERT INTO pending_calculations (mol_hash, model_name, frag_indices, use_cp)
                VALUES (hash, model, indices, True);

            INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
                VALUES (hash, model, indices, '{}', '{}', 'pending', '{}', False);

            INSERT INTO pending_calculations (mol_hash, model_name, frag_indices, use_cp)
                VALUES (hash, model, indices, False);
          ELSE
            INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
                VALUES (hash, model, indices, '{}', '{}', 'pending', '{}', True);

            INSERT INTO pending_calculations (mol_hash, model_name, frag_indices, use_cp)
                VALUES (hash, model, indices, True);
          END IF;

        ELSE

          INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
              VALUES (hash, model, indices, '{}', '{}', 'pending', '{}', False);

          INSERT INTO pending_calculations (mol_hash, model_name, frag_indices, use_cp)
              VALUES (hash, model, indices, False);

        END IF;

      END LOOP;

      INSERT INTO tags VALUES (hash, model, '{}');

      new_config :=  True;
    ELSE
      new_config := False;
    END IF;

    FOREACH tag_name IN ARRAY tags LOOP
      IF NOT EXISTS(SELECT mol_hash FROM tags WHERE mol_hash=hash AND model_name=model AND tag_name=ANY(tag_names)) THEN
        UPDATE tags SET tag_names=(tag_name || tag_names) WHERE mol_hash=hash AND model_name=model;
      END IF;
    END LOOP;

    IF optimized = True THEN
      IF NOT EXISTS(SELECT mol_hash FROM optimized_geometries WHERE mol_name=name AND mol_hash=hash AND model_name=model) THEN
        INSERT INTO optimized_geometries VALUES (name, hash, model);
      END IF;
    END IF;

    RETURN new_config;
  END;

$$;

create function get_1b_training_set(molecule_name character varying, model character varying, input_tags character varying[], batch_offset integer, batch_size integer) returns TABLE(coords double precision[], energy double precision)
	security definer
	SET search_path=public, pg_temp
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
        tag_name varchar;
      BEGIN

        -- Check to make sure the current user has read priveleges on all training sets before allowing them to get
        -- the training set.

        FOREACH tag_name IN ARRAY input_tags
        LOOP
          IF NOT training_set_exists(tag_name)
          THEN
            -- If the training set does not exist, Error
            raise EXCEPTION 'Training set %% does not exist', tag_name;
          ELSIF NOT has_read_privilege(tag_name)
          THEN
            -- If training set does exist and current user doesn't have read privileges, Error.
            raise EXCEPTION 'User %% does not have read privileges on training set %%', session_user, tag_name;
          END IF;
        END LOOP;

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

          IF NOT mol_tags ISNULL THEN
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
          end if;

        END LOOP;
      END;

$$;

create function get_2b_training_set(molecule_name character varying, monomer1_name character varying, monomer2_name character varying, model character varying, input_tags character varying[], batch_offset integer, batch_size integer) returns TABLE(coords double precision[], binding_energy double precision, interaction_energy double precision, monomer1_deformation_energy double precision, monomer2_deformation_energy double precision)
	security definer
	SET search_path=public, pg_temp
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
        monomer1_cp_energies FLOAT[];
        monomer2_cp_energies FLOAT[];
        coordinates FLOAT[];
        tag_name varchar;
      BEGIN

        -- Check to make sure the current user has read priveleges on all training sets before allowing them to get
        -- the training set.

        FOREACH tag_name IN ARRAY input_tags
        LOOP
          IF NOT training_set_exists(tag_name)
          THEN
            -- If the training set does not exist, Error
            raise EXCEPTION 'Training set %% does not exist', tag_name;
          ELSIF NOT has_read_privilege(tag_name)
          THEN
            -- If training set does exist and current user doesn't have read privileges, Error.
            raise EXCEPTION 'User %% does not have read privileges on training set %%', session_user, tag_name;
          END IF;
        END LOOP;

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
            SELECT energies, status FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0, 1}' AND use_cp = False
                INTO dimer_energies, dimer_status;

            IF dimer_status = 'complete' THEN
              SELECT atom_coordinates  FROM molecule_list WHERE mol_hash = hash
                INTO coords;

              SELECT energies FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0}' AND use_cp = False
                INTO monomer1_energies;
              SELECT energies FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{1}' AND use_cp = False
                INTO monomer2_energies;

              IF substring(model, char_length(model) - 3, 4) = 'True' THEN
                SELECT energies FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0}' AND use_cp = True
                  INTO monomer1_cp_energies;
                SELECT energies FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{1}' AND use_cp = True
                  INTO monomer2_cp_energies;

                interaction_energy := dimer_energies[1] - monomer1_cp_energies[1] - monomer2_cp_energies[1];
              ELSE
                interaction_energy := dimer_energies[1] - monomer1_energies[1] - monomer2_energies[1];
              END IF;

              monomer1_deformation_energy := monomer1_energies[1] - optimized_monomer1_energy;
              monomer2_deformation_energy := monomer2_energies[1] - optimized_monomer2_energy;
              binding_energy := interaction_energy + monomer1_deformation_energy + monomer2_deformation_energy;

              RETURN NEXT;
            END IF;
          END IF;

        END LOOP;
      END;

$$;

create function get_failed_configs(molecule_name character varying, model character varying, input_tags character varying[], batch_offset integer, batch_size integer) returns TABLE(coords double precision[], frags integer[], used_cp boolean)
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
        hash VARCHAR;
        mol_tags VARCHAR[];
        mol_tag VARCHAR;
        valid_tags BOOL;
        mol_properties molecule_properties;
        tag_name varchar;
      BEGIN

        -- Check to make sure the current user has read priveleges on all training sets before allowing them to get
        -- the failed configs.

        FOREACH tag_name IN ARRAY input_tags
        LOOP
          IF NOT training_set_exists(tag_name)
          THEN
            -- If the training set does not exist, Error
            raise EXCEPTION 'Training set %% does not exist', tag_name;
          ELSIF NOT has_read_privilege(tag_name)
          THEN
            -- If training set does exist and current user doesn't have read privileges, Error.
            raise EXCEPTION 'User %% does not have read privileges on training set %%', session_user, tag_name;
          END IF;
        END LOOP;

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
              IF mol_properties.status = 'failed' THEN
                SELECT atom_coordinates FROM molecule_list WHERE mol_hash = hash
                  INTO coords;
                frags := mol_properties.frag_indices;
                used_cp = mol_properties.use_cp;
                RETURN NEXT;
              END IF;
          END IF;

        END LOOP;
      END;

$$;

create function get_user_id() returns integer
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  username varchar;
      BEGIN
        SELECT * FROM session_user
            INTO username;

        RETURN get_user_id(username);
      END;

$$;

create function import_calculation(hash character varying, name character varying, fragments fragment[], counts integer[], coordinates double precision[], method character varying, basis character varying, cp boolean, tags character varying[], optimized boolean, nmer_energies double precision[]) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  model varchar;
  indices INTEGER[];
  all_frags INTEGER[];
  i INTEGER;
  x integer;
  z integer;
  new_config BOOLEAN;
  tag_name VARCHAR;
  energy_index INTEGER;
  BEGIN
    model := concat(method, '/', basis, '/');
    IF cp THEN
      model := concat(model, 'True');
    ELSE
      model := concat(model, 'False');
    end if;

    -- Check to Make sure the current user has write priveleges on all training sets specified in the tags.

    FOREACH tag_name IN ARRAY tags
    LOOP
      IF NOT training_set_exists(tag_name)
      THEN
        -- If the training set does not exist, create it giving current user all privileges
        INSERT INTO training_sets (admins, tag_name, read_users, write_users)
            VALUES (ARRAY[get_user_id()], tag_name, ARRAY[get_user_id()], ARRAY[get_user_id()]);
      ELSIF NOT has_write_privilege(tag_name)
      THEN
        -- If training set does exist and current user doesn't have write privileges, Error.
        raise EXCEPTION 'User %% does not have write privileges on training set %%', session_user, tag_name;
      END IF;
    END LOOP;

    IF NOT EXISTS(SELECT mol_hash FROM molecule_properties WHERE mol_hash=hash AND model_name=model) THEN
      PERFORM add_molecule(hash, name, fragments, counts, coordinates);
      PERFORM add_model_info(method, basis, cp);


      i = 0;
      FOREACH x IN ARRAY counts LOOP
        FOR z IN 1..x LOOP
          all_frags := all_frags || i;
          i := i + 1;
        END LOOP;
      END LOOP;

      energy_index := 1;

      FOR indices IN SELECT perm FROM combinations(all_frags) ORDER BY l, perm LOOP

        IF cp = True
        THEN

          IF array_length(indices, 1) = array_length(all_frags, 1)
          THEN

            INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
                VALUES (hash, model, indices, ARRAY[nmer_energies[energy_index]], '{}', 'complete', '{}', False);
            energy_index := energy_index + 1;

          ELSEIF array_length(indices, 1) = 1
          THEN

            INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
                VALUES (hash, model, indices, ARRAY[nmer_energies[energy_index]], '{}', 'complete', '{}', True);
            energy_index := energy_index + 1;

            INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
                VALUES (hash, model, indices, ARRAY[nmer_energies[energy_index]], '{}', 'complete', '{}', False);
            energy_index := energy_index + 1;

          ELSE

            INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
                VALUES (hash, model, indices, ARRAY[nmer_energies[energy_index]], '{}', 'complete', '{}', True);
            energy_index := energy_index + 1;

          END IF;

        ELSE

          INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
              VALUES (hash, model, indices, ARRAY[nmer_energies[energy_index]], '{}', 'complete', '{}', False);
            energy_index := energy_index + 1;

        END IF;

      END LOOP;

      INSERT INTO tags VALUES (hash, model, '{}');

      new_config :=  True;
    ELSE
      new_config := False;
    END IF;

    FOREACH tag_name IN ARRAY tags LOOP
      IF NOT EXISTS(SELECT mol_hash FROM tags WHERE mol_hash=hash AND model_name=model AND tag_name=ANY(tag_names)) THEN
        UPDATE tags SET tag_names=(tag_name || tag_names) WHERE mol_hash=hash AND model_name=model;
      END IF;
    END LOOP;

    IF optimized = True THEN
      IF NOT EXISTS(SELECT mol_hash FROM optimized_geometries WHERE mol_name=name AND mol_hash=hash AND model_name=model) THEN
        INSERT INTO optimized_geometries VALUES (name, hash, model);
      END IF;
    END IF;

    RETURN new_config;
  END;

$$;

create function grant_admin_privilege(input_tag character varying, username character varying) returns void
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    cur_admins int[];
      BEGIN

        -- Check if current user has admin privileges before we allow them to grant admin privileges.
        IF not has_admin_privilege(input_tag)
        THEN
          -- If current user does not have admin privileges, error
          raise EXCEPTION 'User %% does not have admin privileges on training set %%', session_user, input_tag;
        END IF;

      SELECT admins FROM training_sets WHERE tag_name = input_tag
          INTO cur_admins;

      IF NOT get_user_id(username) = ANY(cur_admins)
      THEN
        UPDATE training_sets SET admins = admins || get_user_id(username) WHERE tag_name = input_tag;
      END IF;
    END;

$$;

create function grant_read_privilege(input_tag character varying, username character varying) returns void
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    cur_readers int[];
      BEGIN

        -- Check if current user has admin privileges before we allow them to grant read privileges.
        IF not has_admin_privilege(input_tag)
        THEN
          -- If current user does not have admin privileges, error
          raise EXCEPTION 'User %% does not have admin privileges on training set %%', session_user, input_tag;
        END IF;

      SELECT read_users FROM training_sets WHERE tag_name = input_tag
          INTO cur_readers;

      IF NOT get_user_id(username) = ANY(cur_readers)
      THEN
        UPDATE training_sets SET read_users = read_users || get_user_id(username) WHERE tag_name = input_tag;
      END IF;
    END;

$$;

create function grant_write_privilege(input_tag character varying, username character varying) returns void
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    cur_writers int[];
      BEGIN

        -- Check if current user has admin privileges before we allow them to grant write privileges.
        IF not has_admin_privilege(input_tag)
        THEN
          -- If current user does not have admin privileges, error
          raise EXCEPTION 'User %% does not have admin privileges on training set %%', session_user, input_tag;
        END IF;

      SELECT write_users FROM training_sets WHERE tag_name = input_tag
          INTO cur_writers;

      IF NOT get_user_id(username) = ANY(cur_writers)
      THEN
        UPDATE training_sets SET write_users = write_users || get_user_id(username) WHERE tag_name = input_tag;
      END IF;
    END;

$$;

create function revoke_admin_privilege(input_tag character varying, username character varying) returns void
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    cur_admins int[];
      BEGIN

        -- Check if current user has admin privileges before we allow them to revoke admin privileges.
        IF not has_admin_privilege(input_tag)
        THEN
          -- If current user does not have admin privileges, error
          raise EXCEPTION 'User %% does not have admin privileges on training set %%', session_user, input_tag;
        END IF;

      SELECT admins FROM training_sets WHERE tag_name = input_tag
          INTO cur_admins;

      -- Make sure the the user is currently in the list of admins
      IF NOT get_user_id(username) = ANY(cur_admins)
      THEN
        -- If username is not in the list of admins, error
        raise EXCEPTION 'User %% is already not an admin of training set %%', username, input_tag;
      END IF;

      IF array_length(cur_admins, 1) = 1
      THEN
        -- If the user is trying to remove the last admin, error.
        raise EXCEPTION 'User %% is the last admin of this training set and many not be removed.', username;

      END IF;

      UPDATE training_sets SET admins = array_remove(admins, get_user_id(username)) WHERE tag_name = input_tag;
    END;

$$;

create function revoke_read_privilege(input_tag character varying, username character varying) returns void
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    cur_readers int[];
      BEGIN

        -- Check if current user has admin privileges before we allow them to revoke read privileges.
        IF not has_admin_privilege(input_tag)
        THEN
          -- If current user does not have admin privileges, error
          raise EXCEPTION 'User %% does not have admin privileges on training set %%', session_user, input_tag;
        END IF;

      SELECT read_users FROM training_sets WHERE tag_name = input_tag
          INTO cur_readers;

      -- Make sure the the user is currently in the list of readers.
      IF NOT get_user_id(username) = ANY(cur_readers)
      THEN
        -- If username is not in the list of readers, error.
        raise EXCEPTION 'User %% is already not a reader of training set %%', username, input_tag;
      END IF;

      UPDATE training_sets SET read_users = array_remove(read_users, get_user_id(username)) WHERE tag_name = input_tag;
    END;

$$;

create function revoke_write_privilege(input_tag character varying, username character varying) returns void
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    cur_writers int[];
      BEGIN

        -- Check if current user has admin privileges before we allow them to revoke write privileges.
        IF not has_admin_privilege(input_tag)
        THEN
          -- If current user does not have admin privileges, error
          raise EXCEPTION 'User %% does not have admin privileges on training set %%', session_user, input_tag;
        END IF;

      SELECT write_users FROM training_sets WHERE tag_name = input_tag
          INTO cur_writers;

      -- Make sure the the user is currently in the list of writers.
      IF NOT get_user_id(username) = ANY(cur_writers)
      THEN
        -- If username is not in the list of writers, error.
        raise EXCEPTION 'User %% is already not a writer of training set %%', username, input_tag;
      END IF;

      UPDATE training_sets SET write_users = array_remove(write_users, get_user_id(username)) WHERE tag_name = input_tag;
    END;

$$;

create function delete_calculation(hash character varying, name character varying, method character varying, basis character varying, cp boolean, tags character varying[], delete_complete_calculations boolean) returns boolean
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  model varchar;
  tag_name VARCHAR;
  can_delete BOOLEAN;
  stat VARCHAR;
  BEGIN
    model := concat(method, '/', basis, '/');
    IF cp THEN
      model := concat(model, 'True');
    ELSE
      model := concat(model, 'False');
    end if;

    -- Check to Make sure the current user has write priveleges on all training sets for this geometry

    FOREACH tag_name IN ARRAY tags
    LOOP
      IF NOT training_set_exists(tag_name)
      THEN
        -- If the training set does not exist, Error
        raise EXCEPTION 'Training set %% does not exist', tag_name;
      ELSIF NOT has_write_privilege(tag_name)
      THEN
        -- If training set does exist and current user doesn't have write privileges, Error.
        raise EXCEPTION 'User %% does not have write privileges on training set %%', session_user, tag_name;
      END IF;
    END LOOP;

    IF EXISTS(SELECT mol_hash FROM tags WHERE mol_hash=hash AND model_name=model AND delete_calculation.tags && tags.tag_names) THEN

      FOREACH tag_name in ARRAY tags LOOP
        UPDATE tags SET tag_names=array_remove(tag_names, tag_name) WHERE mol_hash=hash AND model_name=model;
      END LOOP;

      IF EXISTS(SELECT tag_names FROM tags WHERE mol_hash=hash AND model_name=model AND tag_names = '{}') THEN

        DELETE FROM optimized_geometries WHERE mol_hash = hash AND mol_name = name AND model_name = model;

        can_delete := True;

        if not delete_complete_calculations THEN

          for stat in SELECT status FROM molecule_properties WHERE mol_hash=hash AND model_name=model LOOP
            if stat != 'pending' THEN
              can_delete = False;
            end if;
          end loop;

        end if;

        if can_delete THEN
          DELETE FROM tags WHERE mol_hash=hash AND model_name=model;
          DELETE FROM molecule_properties WHERE mol_hash=hash AND model_name=model;
          DELETE FROM pending_calculations WHERE mol_hash=hash AND model_name=model;
          PERFORM delete_molecule(hash, name);
          PERFORM delete_model_info(method, basis, cp);
        end if;


      END IF;


      RETURN TRUE;

    end if;

    RETURN FALSE;
  END;

$$;

create function get_pending_molecule_name(input_tags character varying[]) returns character varying
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
  hash varchar;
  n varchar;
  c integer;
  BEGIN
    SELECT pending_calculations.mol_hash
        FROM pending_calculations
        INNER JOIN tags
        ON pending_calculations.mol_hash = tags.mol_hash AND pending_calculations.model_name = tags.model_name
        WHERE tags.tag_names && input_tags LIMIT 1
    INTO hash;

    IF hash ISNULL
    THEN
      RETURN '';
    END IF;

    SELECT mol_name
        FROM molecule_list
        WHERE mol_hash = hash
    INTO n;

    RETURN n;
  END;

$$;

create function count_training_set_size(molecule_name character varying, model character varying, input_tags character varying[]) returns integer
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    count integer;
    tag_name varchar;
  BEGIN


    -- Check to make sure the current user has read priveleges on all training sets before allowing them to count
    -- the size of the training set.

    FOREACH tag_name IN ARRAY input_tags
    LOOP
      IF NOT training_set_exists(tag_name)
      THEN
        -- If the training set does not exist, Error
        raise EXCEPTION 'Training set %% does not exist', tag_name;
      ELSIF NOT has_read_privilege(tag_name)
      THEN
        -- If training set does exist and current user doesn't have read privileges, Error.
        raise EXCEPTION 'User %% does not have read privileges on training set %%', session_user, tag_name;
      END IF;
    END LOOP;

    SELECT COUNT(*) FROM molecule_list INNER JOIN tags
            ON molecule_list.mol_hash = tags.mol_hash
            WHERE molecule_list.mol_name = molecule_name AND tags.model_name = model AND tags.tag_names && input_tags
            AND 'complete'=ALL(SELECT status FROM molecule_properties WHERE molecule_properties.mol_hash = molecule_list.mol_hash AND molecule_properties.model_name
            = model)
      into count;
    RETURN count;
  END;

$$;

create function export_calculations(molecule_name character varying, model character varying, input_tags character varying[], batch_offset integer, batch_size integer) returns TABLE(coords double precision[], energies double precision[])
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
        hash VARCHAR;
        mol_tags VARCHAR[];
        mol_tag VARCHAR;
        valid_tags BOOL;
        mol_properties molecule_properties;
        tag_name varchar;
      BEGIN

        -- Check to make sure the current user has read priveleges on all training sets before allowing them to export
        -- the configurations in the training set.

        FOREACH tag_name IN ARRAY input_tags
        LOOP
          IF NOT training_set_exists(tag_name)
          THEN
            -- If the training set does not exist, Error
            raise EXCEPTION 'Training set %% does not exist', tag_name;
          ELSIF NOT has_read_privilege(tag_name)
          THEN
            -- If training set does exist and current user doesn't have read privileges, Error.
            raise EXCEPTION 'User %% does not have read privileges on training set %%', session_user, tag_name;
          END IF;
        END LOOP;

        <<main_loop>>
        FOR hash IN SELECT mol_hash FROM molecule_list WHERE mol_name = molecule_name
            ORDER BY mol_hash
            OFFSET batch_offset LIMIT batch_size
        LOOP
          valid_tags := FALSE;


          SELECT tag_names FROM tags WHERE mol_hash = hash AND model_name = model LIMIT 1
            INTO mol_tags;

          IF NOT mol_tags ISNULL THEN
            FOREACH mol_tag IN ARRAY mol_tags
            LOOP


              IF (SELECT mol_tag = ANY(input_tags)) THEN
                valid_tags := TRUE;
              END IF;

            END LOOP;

            IF (SELECT valid_tags) THEN
              energies = '{}';
              for mol_properties IN SELECT * FROM molecule_properties
                  WHERE mol_hash = hash AND model_name = model ORDER BY array_length(frag_indices, 1) ASC, frag_indices ASC, use_cp DESC
              LOOP

                IF mol_properties.status != 'complete' THEN
                  CONTINUE main_loop;
                end if;

                energies = energies || mol_properties.energies[1];

              end loop;

              SELECT atom_coordinates  FROM molecule_list WHERE mol_hash = hash
                    INTO coords;

              RETURN NEXT;
            END IF;
          end if;

        END LOOP;
      END;

$$;

create function count_pending_calculations(input_tags character varying[]) returns integer
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    count integer;
  BEGIN

    SELECT COUNT(*) FROM pending_calculations INNER JOIN tags
            ON pending_calculations.mol_hash = tags.mol_hash
            WHERE tags.tag_names && input_tags
      into count;
    RETURN count;
  END;

$$;

create function get_training_set(molecule_name character varying, monomer_names character varying[], model character varying, input_tags character varying[], batch_offset integer, batch_size integer) returns TABLE(coords double precision[], binding_energy double precision, nb_energy double precision, deformation_energies double precision[])
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
        monomer_name varchar;
        optimized_energies FLOAT[];
        optimized_index INT;
        hash VARCHAR;
        energy FLOAT;
        stat VARCHAR;
        deformation_index INT;
        cp BOOLEAN;
        n_mer INT;
        tag_name varchar;
      BEGIN

        -- Check to make sure the current user has read priveleges on all training sets before allowing them to get
        -- the failed configs.

        FOREACH tag_name IN ARRAY input_tags
        LOOP
          IF NOT training_set_exists(tag_name)
          THEN
            -- If the training set does not exist, Error
            raise EXCEPTION 'Training set %% does not exist', tag_name;
          ELSIF NOT has_read_privilege(tag_name)
          THEN
            -- If training set does exist and current user doesn't have read privileges, Error.
            raise EXCEPTION 'User %% does not have read privileges on training set %%', session_user, tag_name;
          END IF;
        END LOOP;

        -- monomer_names must be passed in in standard order!
        -- First, get each of the optimized energies

        optimized_energies := '{}';
        optimized_index := 1;

        FOREACH monomer_name IN ARRAY monomer_names
        LOOP

          optimized_energies = optimized_energies || NULL;


          FOR hash IN SELECT optimized_geometries.mol_hash
              FROM optimized_geometries INNER JOIN tags
              ON optimized_geometries.mol_hash = tags.mol_hash AND optimized_geometries.model_name = tags.model_name
              WHERE optimized_geometries.mol_name = monomer_name AND optimized_geometries.model_name = model AND tags.tag_names && input_tags
          LOOP

            SELECT energies[1], status FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0}'
              INTO energy, stat;

            IF stat = 'complete'
            THEN

              IF (SELECT optimized_energies[optimized_index] ISNULL)
              THEN

                optimized_energies[optimized_index] := energy;

              ELSE

                IF (energy < optimized_energies[optimized_index])
                THEN

                  optimized_energies[optimized_index] := energy;

                end if;

                RAISE WARNING 'Multiple optimized geometries for %% in database.', monomer_name;

              END IF;

            ELSE

              RAISE EXCEPTION 'Optimized energy for %% uncalculated in database.', monomer_name;

            END IF;



          END LOOP;

          IF optimized_energies[optimized_index] ISNULL
          THEN

            RAISE EXCEPTION 'No optimized energy in database for %%.', monomer_name;

          END IF;


          optimized_index = optimized_index + 1;


        END LOOP;

        -- Does the model use cp?

        IF substring(model, char_length(model) - 3, 4) = 'True'
        THEN
          cp = True;
        ELSE
          cp = False;
        END IF;

        -- LOOP over each element of the training set.

        <<TrainingSetLoop>>
        FOR hash IN SELECT molecule_list.mol_hash, tags.mol_hash
            FROM molecule_list INNER JOIN tags
            ON molecule_list.mol_hash = tags.mol_hash
            WHERE molecule_list.mol_name = molecule_name AND tags.model_name = model AND tags.tag_names && input_tags
            AND 'complete'=ALL(SELECT status FROM molecule_properties WHERE molecule_properties.mol_hash = molecule_list.mol_hash AND molecule_properties.model_name
            = model)
            ORDER BY molecule_list.mol_hash
            OFFSET batch_offset LIMIT batch_size
        LOOP
          -- Get the deformation energy for each monomer
          deformation_energies := '{}';
          deformation_index := 1;
          FOR energy, stat in SELECT energies[1], status
              FROM molecule_properties
              WHERE mol_hash = hash AND model_name = model AND array_length(frag_indices, 1) = 1 AND use_cp = False
              ORDER BY frag_indices ASC
          LOOP
            IF stat != 'complete'
            THEN
              raise EXCEPTION 'NEVER';
              CONTINUE TrainingSetLoop;
            END IF;
            deformation_energies := deformation_energies || energy - optimized_energies[deformation_index];
            deformation_index := deformation_index + 1;
          END LOOP;

          -- Calculate the nb_energy energy.

          IF array_length(monomer_names, 1) = 1
          THEN
            -- For monomers, the nb_energy energy and binding energy are both the monomer deformation energy.
            nb_energy := deformation_energies[1];
            binding_energy := nb_energy;
          ELSE
            -- For non-monomers, calculate the nb_energy energy normally.
            SELECT energies[1], status
                FROM molecule_properties
                WHERE mol_hash = hash AND model_name = model AND array_length(frag_indices, 1) = array_length(monomer_names, 1)
            INTO nb_energy, stat;

            IF stat != 'complete'
            THEN
              raise EXCEPTION 'NEVER';
              CONTINUE TrainingSetLoop;
            END IF;

            FOR energy, stat, n_mer IN SELECT energies[1], status, array_length(frag_indices, 1)
                FROM molecule_properties
                WHERE mol_hash = hash AND model_name = model AND array_length(frag_indices, 1) < array_length(monomer_names, 1) AND use_cp = cp
            LOOP

              IF stat != 'complete'
              THEN
              raise EXCEPTION 'NEVER';
                CONTINUE TrainingSetLoop;
              END IF;

              IF (array_length(monomer_names, 1) - n_mer) %% 2 = 1
              THEN
                nb_energy := nb_energy - energy;
              ELSE
                nb_energy := nb_energy + energy;
              END IF;

            END LOOP;

            -- For non-monomers, calculate the binding energy normally.

            binding_energy := nb_energy;
            FOREACH energy IN ARRAY deformation_energies
            LOOP
              binding_energy := binding_energy + energy;
            END LOOP;
          END IF;


          SELECT atom_coordinates
              FROM molecule_list
              WHERE mol_hash = hash
          INTO coords;

          RETURN NEXT;

        END LOOP;

      END;
$$;

create function count_dispatched_calculations() returns integer
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    count integer;
  BEGIN

    SELECT COUNT(*) FROM molecule_properties INNER JOIN tags
            ON molecule_properties.mol_hash = tags.mol_hash
            WHERE molecule_properties.status = 'dispatched'
      into count;
    RETURN count;
  END;

$$;

create function set_properties(hash character varying, model character varying, use_cp boolean, indices integer[], result boolean, energy double precision, log_txt character varying, overwrite boolean) returns void
	security definer
	SET search_path=public, pg_temp
	language plpgsql
as $$
DECLARE
    id INTEGER;
    cur_status VARCHAR;
  BEGIN

    SELECT status FROM molecule_properties WHERE  mol_hash=hash AND model_name=model AND frag_indices=indices AND molecule_properties.use_cp = set_properties.use_cp
      INTO cur_status;

    IF cur_status = 'complete' THEN
      IF overwrite THEN
        UPDATE molecule_properties SET energies = '{}', status='dispatched' WHERE mol_hash=hash AND model_name=model AND frag_indices=indices AND molecule_properties.use_cp = set_properties.use_cp;
        cur_status = 'dispatched';
      END IF;
    END IF;

    IF NOT cur_status = 'dispatched' THEN
      IF cur_status = 'complete' THEN
        RETURN;
      END IF;
      RAISE EXCEPTION 'Trying to set energy of calculation that does not have status = "dispatched" or status = "complete"';
    END IF;

    IF result THEN
      UPDATE molecule_properties SET energies = energy || energies, status='complete' WHERE mol_hash=hash AND model_name=model AND frag_indices=indices AND molecule_properties.use_cp = set_properties.use_cp RETURNING most_recent_log_id
        INTO id;
      UPDATE log_files SET end_time=clock_timestamp(), log_text=log_txt WHERE log_id=id;
    ELSE
      UPDATE molecule_properties SET energies = energies || energy, status='failed' WHERE mol_hash=hash AND model_name=model AND frag_indices=indices AND molecule_properties.use_cp = set_properties.use_cp RETURNING most_recent_log_id
        INTO id;
      UPDATE log_files SET end_time=clock_timestamp(), log_text=log_txt WHERE log_id=id;
    END IF;


  END;

$$;

