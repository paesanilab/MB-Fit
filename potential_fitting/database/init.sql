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

create function set_properties(hash character varying, model character varying, indices integer[], result boolean, energy double precision, log_txt character varying) returns void
	language plpgsql
as $$
DECLARE
    id INTEGER;
    cur_status VARCHAR;
  BEGIN

    SELECT status FROM molecule_properties WHERE  mol_hash=hash AND model_name=model AND frag_indices=indices
      INTO cur_status;
    IF NOT cur_status = 'dispatched' THEN
      RAISE EXCEPTION 'Trying to set energy of calculation that does not have status = "dispatched"';
    END IF;

    IF result THEN
      UPDATE molecule_properties SET energies = energies || energy, status='complete' WHERE mol_hash=hash AND model_name=model AND frag_indices=indices RETURNING most_recent_log_id
        INTO id;
      UPDATE log_files SET end_time=clock_timestamp(), log_text=log_txt WHERE log_id=id;
    ELSE
      UPDATE molecule_properties SET energies = energies || energy, status='failed' WHERE mol_hash=hash AND model_name=model AND frag_indices=indices RETURNING most_recent_log_id
        INTO id;
      UPDATE log_files SET end_time=clock_timestamp(), log_text=log_txt WHERE log_id=id;
    END IF;


  END;

$$;

alter function set_properties(varchar, varchar, integer[], boolean, double precision, varchar) owner to ebullvul;

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

alter function count_entries(varchar) owner to potential_fitting;