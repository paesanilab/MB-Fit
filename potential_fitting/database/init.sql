create schema public;

comment on schema public is 'standard public schema';

alter schema public owner to ebullvul;

create type empty_mol as
(
	mol_name varchar,
	frag_names character varying[],
	frag_counts integer[],
	frag_charges integer[],
	frag_spins integer[],
	frag_smiles character varying[]
);

alter type empty_mol owner to ebullvul;

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

alter type fragment owner to ebullvul;

create type status_enum as enum ('pending', 'dispatched', 'complete', 'failed');

alter type status_enum owner to ebullvul;

create table molecule_info
(
	name varchar not null
		constraint "molecule_Name_index"
			primary key
);

alter table molecule_info owner to ebullvul;

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

alter table molecule_list owner to ebullvul;

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

alter table atom_info owner to ebullvul;

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

alter table fragment_info owner to ebullvul;

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

alter table molecule_contents owner to ebullvul;

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

alter table fragment_contents owner to ebullvul;

create index fragment_contents_frag_name_index
	on fragment_contents (frag_name);

create table model_info
(
	name varchar not null
		constraint models_info_name_key
			primary key
);

alter table model_info owner to ebullvul;

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

alter table log_files owner to ebullvul;

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

alter table molecule_properties owner to ebullvul;

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

alter table tags owner to ebullvul;

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

alter table optimized_geometries owner to ebullvul;

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

alter table pending_calculations owner to ebullvul;

create index pending_calculations_mol_hash_model_name_frag_indices_use_cp_in
	on pending_calculations (mol_hash, model_name, frag_indices, use_cp);

create function get_empty_molecule(molecule_name character varying) returns empty_mol
	language plpgsql
as $$
DECLARE
    mol empty_mol;
    fragment_desc RECORD;
  BEGIN
  mol.mol_name := molecule_name;
  FOR fragment_desc IN SELECT frag_name, count, charge, spin, SMILE
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
      mol.frag_SMILES := mol.frag_SMILES || fragment_desc.SMILE;

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

alter function get_1b_training_set(varchar, varchar, character varying[], integer, integer) owner to ebullvul;

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
        monomer1_cp_energies FLOAT[];
        monomer2_cp_energies FLOAT[];
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
            SELECT energies, status FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0, 1}' AND use_cp = False
                INTO dimer_energies, dimer_status;

            IF dimer_status = 'complete' THEN
              SELECT atom_coordinates  FROM molecule_list WHERE mol_hash = hash
                INTO coords;

              SELECT energies FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{0}' AND use_cp = False
                INTO monomer1_energies;
              SELECT energies FROM molecule_properties WHERE mol_hash = hash AND model_name = model AND frag_indices = '{1}' AND use_cp = False
                INTO monomer2_energies;

              IF substring(model, char_length(model) - 4, 4) = 'True' THEN
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

alter function get_2b_training_set(varchar, varchar, varchar, varchar, character varying[], integer, integer) owner to ebullvul;

create function count_entries(molecule_name character varying) returns integer
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

alter function count_entries(varchar) owner to ebullvul;

create function add_atom_info(atomic_symbol character varying) returns boolean
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

alter function add_atom_info(varchar) owner to ebullvul;

create function add_molecule_info(name character varying, fragments fragment[], counts integer[]) returns boolean
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

alter function add_molecule_info(varchar, fragment[], integer[]) owner to ebullvul;

create function add_model_info(method character varying, basis character varying, cp boolean) returns boolean
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

alter function add_model_info(varchar, varchar, boolean) owner to ebullvul;

create function add_molecule(hash character varying, name character varying, fragments fragment[], counts integer[], coordinates double precision[]) returns boolean
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

alter function add_molecule(varchar, varchar, fragment[], integer[], double precision[]) owner to ebullvul;

create function set_properties(hash character varying, model character varying, use_cp boolean, indices integer[], result boolean, energy double precision, log_txt character varying) returns void
	language plpgsql
as $$
DECLARE
    id INTEGER;
    cur_status VARCHAR;
  BEGIN

    SELECT status FROM molecule_properties WHERE  mol_hash=hash AND model_name=model AND frag_indices=indices AND molecule_properties.use_cp = set_properties.use_cp
      INTO cur_status;
    IF NOT cur_status = 'dispatched' THEN
      RAISE EXCEPTION 'Trying to set energy of calculation that does not have status = "dispatched"';
    END IF;

    IF result THEN
      UPDATE molecule_properties SET energies = energies || energy, status='complete' WHERE mol_hash=hash AND model_name=model AND frag_indices=indices AND molecule_properties.use_cp = set_properties.use_cp RETURNING most_recent_log_id
        INTO id;
      UPDATE log_files SET end_time=clock_timestamp(), log_text=log_txt WHERE log_id=id;
    ELSE
      UPDATE molecule_properties SET energies = energies || energy, status='failed' WHERE mol_hash=hash AND model_name=model AND frag_indices=indices AND molecule_properties.use_cp = set_properties.use_cp RETURNING most_recent_log_id
        INTO id;
      UPDATE log_files SET end_time=clock_timestamp(), log_text=log_txt WHERE log_id=id;
    END IF;


  END;

$$;

alter function set_properties(varchar, varchar, boolean, integer[], boolean, double precision, varchar) owner to ebullvul;

create function get_pending_calculations(molecule_name character varying, input_client_name character varying, input_tags character varying[], batch_size integer) returns TABLE(coords double precision[], model character varying, indices integer[], use_cp boolean)
	language plpgsql
as $$
DECLARE
    molecule_hash VARCHAR;
    id INTEGER;

  BEGIN

    FOR molecule_hash, model, indices, get_pending_calculations.use_cp IN SELECT pending_calculations.mol_hash, pending_calculations.model_name, pending_calculations.frag_indices, pending_calculations.use_cp FROM
      pending_calculations INNER JOIN molecule_list ON pending_calculations.mol_hash = molecule_list.mol_hash WHERE
      molecule_list.mol_name = molecule_name LIMIT batch_size
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

alter function get_pending_calculations(varchar, varchar, character varying[], integer) owner to ebullvul;

create function add_calculation(hash character varying, name character varying, fragments fragment[], counts integer[], coordinates double precision[], method character varying, basis character varying, cp boolean, tags character varying[], optimized boolean) returns boolean
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
        IF cp = True AND array_length(indices, 1) != array_length(all_frags, 1) THEN
          INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
              VALUES (hash, model, indices, '{}', '{}', 'pending', '{}', True);

          INSERT INTO pending_calculations (mol_hash, model_name, frag_indices, use_cp)
              VALUES (hash, model, indices, True);
        END IF;

        INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
              VALUES (hash, model, indices, '{}', '{}', 'pending', '{}', False);

          INSERT INTO pending_calculations (mol_hash, model_name, frag_indices, use_cp)
              VALUES (hash, model, indices, False);

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

alter function add_calculation(varchar, varchar, fragment[], integer[], double precision[], varchar, varchar, boolean, character varying[], boolean) owner to ebullvul;

create function combinations(arr integer[]) returns TABLE(perm integer[], l integer)
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

alter function combinations(integer[]) owner to ebullvul;

create function import_calculation(hash character varying, name character varying, fragments fragment[], counts integer[], coordinates double precision[], method character varying, basis character varying, cp boolean, tags character varying[], optimized boolean, nmer_energies double precision[]) returns boolean
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
        IF cp = True AND array_length(indices, 1) != array_length(all_frags, 1) THEN
          INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
              VALUES (hash, model, indices, ARRAY[nmer_energies[energy_index]], '{}', 'complete', '{}', True);

          energy_index := energy_index + 1;
        END IF;

        INSERT INTO molecule_properties (mol_hash, model_name, frag_indices, energies, atomic_charges, status, past_log_ids, use_cp)
              VALUES (hash, model, indices, ARRAY[nmer_energies[energy_index]], '{}', 'complete', '{}', False);

        energy_index := energy_index + 1;
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

alter function import_calculation(varchar, varchar, fragment[], integer[], double precision[], varchar, varchar, boolean, character varying[], boolean, double precision[]) owner to ebullvul;

create function get_failed_configs(molecule_name character varying, model character varying, input_tags character varying[], batch_offset integer, batch_size integer) returns TABLE(coords double precision[], frags integer[], used_cp boolean)
	language plpgsql
as $$
DECLARE
        hash VARCHAR;
        mol_tags VARCHAR[];
        mol_tag VARCHAR;
        valid_tags BOOL;
        mol_properties molecule_properties;
      BEGIN

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

alter function get_failed_configs(varchar, varchar, character varying[], integer, integer) owner to ebullvul;

create function remove_calculation(hash character varying, name character varying, method character varying, basis character varying, cp boolean, tags character varying[], optimized boolean) returns boolean
	language plpgsql
as $$
DECLARE
  model varchar;
  cur_tags varchar[];
  t varchar;
  mol_frags varchar[];
  frag varchar;
  mol_atoms varchar[];
  atom varchar;
  BEGIN
    model := concat(method, '/', basis, '/');
    IF cp THEN
      model := concat(model, 'True');
    ELSE
      model := concat(model, 'False');
    end if;

    -- still must clear optimized_geometries, pending_calculations, and log_files tables.

    IF EXISTS(SELECT tag_names FROM tags WHERE mol_hash = hash AND model_name = model) THEN
      SELECT tag_names FROM tags WHERE mol_hash = hash AND model_name = model
        INTO cur_tags;
      FOREACH t IN ARRAY tags LOOP
        cur_tags = array_remove(cur_tags, t);
      end loop;
      UPDATE tags SET tag_names = cur_tags WHERE mol_hash = hash AND model_name = model;

      IF array_length(CURTAGS, 1) = 0 THEN
        DELETE FROM tags WHERE mol_hash = hash AND model_name = model;
        DELETE FROM molecule_properties WHERE mol_hash = hash and model_name = model;

        IF NOT EXISTS(SELECT mol_hash FROM molecule_properties WHERE mol_hash = hash) THEN
          DELETE FROM molecule_list WHERE mol_hash = hash;

          IF NOT EXISTS(SELECT mol_hash FROM molecule_list WHERE mol_name = name) THEN
            DELETE FROM molecule_info WHERE molecule_info.name = remove_calculation.name;
            DELETE FROM molecule_contents WHERE molecule_contents.mol_name = remove_calculation.name RETURNING frag_name
              INTO mol_frags;

            FOREACH frag in ARRAY mol_frags LOOP
              IF NOT EXISTS(SELECT frag_name FROM molecule_contents WHERE frag_name = frag) THEN
                DELETE FROM fragment_info WHERE name = frag;
                DELETE FROM fragment_contents WHERE frag_name = frag RETURNING atom_symbol
                  INTO mol_atoms;

                FOREACH atom IN ARRAY mol_atoms LOOP
                  IF NOT EXISTS(SELECT atom_symbol FROM fragment_contents WHERE atom_symbol = atom) THEN
                    DELETE FROM atom_info WHERE atomic_symbol = atom;
                  end if;
                end loop;

              end if;
            end loop;


          end if;

        end if;

        IF NOT EXISTS(SELECT moldel_name FROM molecule_properties WHERE model_name = model) THEN
          DELETE FROM model_info WHERE model_info.name = model;
        end if;

      end if;
    end if;
  END;

$$;

alter function remove_calculation(varchar, varchar, varchar, varchar, boolean, character varying[], boolean) owner to ebullvul;

create function reset_dispatched(ts character varying[]) returns integer
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

alter function reset_dispatched(character varying[]) owner to ebullvul;

create function reset_failed(ts character varying[]) returns integer
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

alter function reset_failed(character varying[]) owner to ebullvul;

create function add_fragment_info(name character varying, charge integer, spin integer, smile character varying, atomic_symbols character varying[], symmetries character varying[], counts integer[]) returns boolean
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

alter function add_fragment_info(varchar, integer, integer, varchar, character varying[], character varying[], integer[]) owner to ebullvul;

create function construct_fragment(name character varying, charge integer, spin integer, smile character varying, atomic_symbols character varying[], symmetries character varying[], counts integer[]) returns fragment
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

alter function construct_fragment(varchar, integer, integer, varchar, character varying[], character varying[], integer[]) owner to ebullvul;

