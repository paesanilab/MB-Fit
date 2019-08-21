# external package imports
import os, sys
from hashlib import sha1

# absolute module imports
from potential_fitting.utils import SettingsReader, files, system
from potential_fitting.exceptions import ConfigMissingSectionError, ConfigMissingPropertyError

# local module imports
from .database import Database


def make_all_jobs(settings_path, database_config_path, client_name, job_dir, *tags, num_jobs=sys.maxsize):
    """
    Makes a Job file for each energy that still needs to be calculated in this Database.

    Args:
        settings_path       - Local path to the ".ini" file with relevent settings.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        client_name         - Name of the client that will perform these jobs
        job_dir             - Local path to the directory to place the job files in.
        tags                - Onlt  make jobs for calculations marked with at least one of these tags.
        num_jobs            - The number of jobs to generate. Unlimted if None.

    Returns:
        None.
    """


    if num_jobs is None:
        num_jobs = sys.maxsize

    counter = 0

    # open the database
    with Database(database_config_path) as database:

        total_pending = database.count_pending_calculations(*tags)
        system.format_print("Making jobs from database into directory {}. {} total jobs pending in database. Making jobs for {} of them.".format(job_dir, total_pending, min(num_jobs, total_pending)), bold=True, color=system.Color.YELLOW)

        for molecule, method, basis, cp, use_cp, frag_indices in database.get_all_calculations(client_name, *tags, calculations_to_do=num_jobs):

            write_job(settings_path, molecule, method, basis, cp, use_cp, frag_indices, job_dir)
            counter += 1
            if counter % 100 == 0:
                system.format_print("Made {} jobs so far.".format(counter), italics=True)

    system.format_print("Completed job generation. {} jobs generated. {} jobs remaining to be created.".format(counter, total_pending - counter), bold=True, color=system.Color.GREEN)


def write_job(settings_path, molecule, method, basis, cp, use_cp, frag_indices, job_dir):
    """
    Makes a Job file for a specific calculation.

    cp is not the same as use_cp. Some models have cp, but should not
    use cp for some of their energies.

    Args:
        settings_path       - Local path to the ".ini" file with relevent settings
        molecule            - The molecule of this calculation.
        method              - Method to use to calculate the energy.
        basis               - Basis to use to calculate the energy.
        cp                  - True if the model has counterpoise correction.
        use_cp              - True if counterpoise correction should be used for this calculation.
        frag_indices        - List of indices of fragments to include in the calculation.
        job_dir             - Local path to the directory to place the job file in.

    Returns:
        None.
    
    """

    # parse settings file
    settings = SettingsReader(settings_path)

    template_dictionary = {
            # TODO
            "whole_molecule": molecule.to_xyz().replace("\n", "\\n"),
            "charges": [frag.get_charge() for frag in molecule.get_fragments()],
            "spins": [frag.get_spin_multiplicity() for frag in molecule.get_fragments()],
            "symmetries": [frag.get_symmetry() for frag in molecule.get_fragments()],
            "SMILES": ",".join([frag.get_SMILE() for frag in molecule.get_fragments()]),
            "atom_counts": [frag.get_num_atoms() for frag in molecule.get_fragments()],
            "names": [frag.get_name() for frag in molecule.get_fragments()],
            "total_atoms": molecule.get_num_atoms(),
            "molecule":     molecule.to_xyz(frag_indices, use_cp).replace("\n", "\\n"),
            "frag_indices": frag_indices,
            "method":       method,
            "basis":        basis,
            "cp":           cp,
            "use_cp":       use_cp,
            "num_threads":  settings.get("psi4", "num_threads"),
            "memory":       settings.get("psi4", "memory"),
            "format":       "{}",
            "total_charge": molecule.get_charge(frag_indices),
            "total_spin": molecule.get_spin_multiplicity(frag_indices)
        }

    hash_string = "\n".join([str(v) for v in template_dictionary.values()])
    job_hash = sha1(hash_string.encode()).hexdigest()

    template_dictionary["job_hash"] = job_hash

    i = 8
    file_path = job_dir + "/job_{}.py".format(job_hash[:i])

    while os.path.exists(file_path):
        i += 1

        file_path = job_dir + "/job_{}.py".format(job_hash[:i])

    files.init_file(file_path)

    with open(file_path, "w") as job_file, open(os.path.dirname(os.path.abspath(__file__)) + "/job_template.py", "r") as job_template:
        job_string = "".join(job_template.readlines())

        job_file.write(job_string.format(**template_dictionary))
