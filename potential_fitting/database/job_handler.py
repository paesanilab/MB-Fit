import sys, os
from glob import glob
from hashlib import sha1

from potential_fitting.calculator import Model
from potential_fitting.utils import SettingsReader, system, files
from potential_fitting.molecule import Molecule
from potential_fitting.exceptions import ConfigMissingPropertyError
from . import Database


class JobHander(object):

    def __init__(self, settings_path):
        """
        Constructs a new JobHander.

        Args:
            settings_path   - Local path to '.ini' settings file with all settings.

        Returns:
            A new JobHandler object.
        """
        self.settings = SettingsReader(settings_path)

    def make_all_jobs(self, database_config_path, client_name, job_dir, *tags, num_jobs=sys.maxsize):
        """
        Makes a Job file for each energy that still needs to be calculated in this Database.

        Args:
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
            system.format_print(
                "Making jobs from database into directory {}. {} total jobs with tags {} pending in database. Making jobs for {} of them.".format(
                    job_dir, total_pending, tags, min(num_jobs, total_pending)), bold=True, color=system.Color.YELLOW)

            for molecule, method, basis, cp, use_cp, frag_indices in database.get_all_calculations(client_name, *tags,
                                                                                                   calculations_to_do=num_jobs):
                model = Model(method, basis, cp)
                self.write_job(molecule, model, use_cp, frag_indices, job_dir)
                counter += 1
                if counter % 100 == 0:
                    system.format_print("Made {} jobs so far.".format(counter), italics=True)

        system.format_print(
            "Completed job generation. {} jobs generated. {} jobs with tags {} remaining to be created.".format(counter,
                                                                                                                total_pending - counter,
                                                                                                                tags),
            bold=True, color=system.Color.GREEN)

    def read_all_jobs(self, database_config_path, job_dir):
        """
        Searches the given directory for completed job directories and enters
        the results into the database.

        Any directory that starts with job_ will be considered a completed job directory.
        After data is entered into the database, these directories will be renamed from
        job_* to job_*_done.

        Args:
            database_config_path - .ini file containing host, port, database, username, and password.
                        Make sure only you have access to this file or your password will be compromised!
            job_dir             - Local path the the directory to search.

        Returns:
            None.
        """

        system.format_print("Reading jobs from directory {} into database.".format(job_dir),
                            bold=True, color=system.Color.YELLOW)

        counter = 0

        calculation_results = []
        for directory in glob(job_dir + "/job_*"):
            if directory.endswith("done"):
                continue
            if not os.path.isdir(directory):
                continue
            calculation_results.append(self.read_job(directory + "/output.ini", directory + "/output.log"))

            if len(calculation_results) > 1000:
                with Database(database_config_path) as db:
                    db.set_properties(calculation_results)
                calculation_results = []

            counter += 1

            if counter % 100 == 0:
                system.format_print("Read {} jobs so far.".format(counter), italics=True)

        with Database(database_config_path) as db:
            db.set_properties(calculation_results)

        system.format_print("Completed reading jobs. Read {} in total".format(counter), bold=True,
                            color=system.Color.GREEN)
        system.format_print("Moving read job directories...", bold=True, color=system.Color.YELLOW)

        for directory in glob(job_dir + "/job_*"):
            if directory.endswith("done"):
                continue
            if not os.path.isdir(directory):
                continue

            i = 1

            job_dir = directory + "_done"

            while os.path.exists(job_dir):
                i += 1

                job_dir = directory + "_{}_done".format(i)

            os.rename(directory, job_dir)

        system.format_print("Moved completed job directories!", bold=True, color=system.Color.GREEN)

    def write_job(self, molecule, model, use_cp, frag_indices, job_dir):
        """
        Makes a Job file for a specific calculation.

        cp is not the same as use_cp. Some models have cp, but should not
        use cp for some of their energies.

        Args:
            molecule            - The molecule of this calculation.
            model               - The model to use for this calculation.
            use_cp              - True if counterpoise correction should be used for this calculation.
            frag_indices        - List of indices of fragments to include in the calculation.
            job_dir             - Local path to the directory to place the job file in.

        Returns:
            None.
        """

        # parse settings file

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
                "method":       model.get_method(),
                "basis":        model.get_basis(),
                "cp":           model.get_cp(),
                "use_cp":       use_cp,
                "num_threads":  self.settings.get("psi4", "num_threads"),
                "memory":       self.settings.get("psi4", "memory"),
                "total_charge": molecule.get_charge(frag_indices),
                "total_spin": molecule.get_spin_multiplicity(frag_indices),
                "format":       "{}"
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

        with open(file_path, "w") as job_file, open(self.get_job_template_path(), "r") as job_template:
            job_string = "".join(job_template.readlines())

            job_file.write(job_string.format(**template_dictionary))

    def read_job(self, job_dat_path, job_log_path):
        """
        Reads a completed job from its output file and log_file and enters the result into a database.

        Args:
            job_dat_path        - Local path to the .ini output file from this job.
            job_log_path        - Local path to the log file from this job.

        Returns:
            None.
        """

        data = SettingsReader(job_dat_path)

        xyz = data.get("molecule", "xyz")
        atom_counts = data.getlist("molecule", "atom_counts", int)
        charges = data.getlist("molecule", "charges", int)
        spins = data.getlist("molecule", "spins", int)
        symmetries = data.getlist("molecule", "symmetries", str)
        SMILES = data.get("molecule", "SMILES", str).split(",")
        names = data.getlist("molecule", "names", str)
        method = data.get("molecule", "method")
        basis = data.get("molecule", "basis")
        cp = data.get("molecule", "cp")
        use_cp = data.get("molecule", "use_cp")
        frag_indices = data.getlist("molecule", "frag_indices", int)

        molecule = Molecule.read_xyz(xyz, atom_counts, names, charges, spins, symmetries, SMILES)

        try:
            energy = data.getfloat("molecule", "energy")
            success = True
        except ConfigMissingPropertyError:
            energy = 0
            success = False

        with open(job_log_path, "r") as log_file:
            log_text = "\n".join(log_file.readlines())

        return molecule, method, basis, cp, use_cp, frag_indices, success, energy, log_text

    def get_job_template_path(self):
        raise NotImplementedError
