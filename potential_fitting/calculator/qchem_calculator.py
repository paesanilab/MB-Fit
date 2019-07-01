# absolute module imports
from potential_fitting.utils import files, system
from potential_fitting.molecule import Molecule
from potential_fitting.exceptions import (LibraryNotAvailableError, LibraryCallError,
        ConfigMissingSectionError, ConfigMissingPropertyError, CommandExecutionError)

from . import Calculator

class QchemCalculator(Calculator):

    def __init__(self, settings_path, logging = True):
        """
        Constructor for a QchemCalculator.

        A QchemCalculator uses Qchem for all of its calculations.

        Args:
            settings_path   - Path to the settings.ini
            logging         - Whether this calculator should output logging messages.
                Default: True

        Returns:
            None
        """
        
        super(QchemCalculator, self).__init__(settings_path, logging)

        if not self.is_installed:
            raise LibraryNotAvailableError("qchem")

    def is_installed(self):
        try:
            system.call("which", "qchem")
            return True
        except CommandExecutionError:
            return False

    def create_input_file(self, file_path, molecule, model, job, fragment_indicies = None):
        """
        Creates an input file for a Qchem calculation in the given file path.

        Args:
            file_path       - The path to the file to write the input into.
            molecule        - The molecule that will be used in the calculation.
            model           - The model that will be used in the calculation. If the model has cp enabled, then fragments
                    not included in fragment_indicies will be included as ghost atoms.
            job             - What job this input file should specify. For example: "sp", "opt", or "freq".
            fragment_indicies - If not None, then only these fragments will be included in the input file.
                Defualt: None

        Returns:
            None
        """

        with open(file_path, "w") as in_file:
            in_file.write("$molecule\n")

            in_file.write("{} {}\n".format(molecule.get_charge(fragment_indicies), molecule.get_spin_multiplicity(fragment_indicies)))

            in_file.write("{}\n".format(molecule.to_xyz(fragment_indicies, model.get_cp())))

            in_file.write("$end\n")

            in_file.write("\n")

            in_file.write("$rem\n")

            in_file.write("jobtype {}\n".format(job))

            in_file.write("method {}\n".format(model.get_method()))
            in_file.write("basis {}\n".format(model.get_basis()))

            in_file.write("GEOM_OPT_MAX_CYCLES 200\n")

            try:
                in_file.write("ecp {}\n".format(self.settings.get("config_generator", "ecp")))
            except (ConfigMissingSectionError, ConfigMissingPropertyError):
                pass

            in_file.write("$end\n")

    def calculate_energy(self, molecule, model, fragment_indicies):
        """
        Calculates the energy of a subset of the fragments of a molecule with a provided model.

        Fragments not included in fragment_indicies will be included as ghost atoms if model's
        cp is enabled.

        Args:
            molecule        - The Molecule to perform the calculation on.
            fragment_indicies - List of the indicies of the fragments to include in the calculation.
            model           - The Model to use for the calculation.

        Returns:
            (calculated energy, path to the log file)
        """

        if self.logging:
            print("Beginning energy calculation of {} with model {}/{} and cp = {}".format(molecule.get_name(), model.get_method(), model.get_basis(), model.get_cp()))

        # file to write qchem input in
        qchem_in_path = files.get_energy_log_path(self.settings.get("files", "log_path"), molecule, model.get_method(), model.get_basis(), model.get_cp(), "in")
        
        self.create_input_file(qchem_in_path, molecule, model, "sp", fragment_indicies)

        # file to write qchem output in
        qchem_out_path = files.get_energy_log_path(self.settings.get("files", "log_path"), molecule, model.get_method(), model.get_basis(), model.get_cp(), "out")

        # get number of threads
        num_threads = self.settings.getint("qchem", "num_threads", 1)

        # perform system call to run qchem
        try:
            system.call("qchem", "-nt", str(num_threads), qchem_in_path, qchem_out_path)
        except CommandExecutionError as e:
            raise LibraryCallError("qchem", "energy calculation", str(e), log_path=qchem_out_path)

        energy = self.find_energy_in_energy_calculation_output_file(qchem_out_path)

        if self.logging:
            print("Successfully completed energy calculation.")

        return energy, qchem_out_path

    def find_energy_in_energy_calculation_output_file(self, qchem_out_path):
        """
        Parses the output file of a Qchem energy calculation and retrieves the energy that
        resulted from the calculation.

        Args:
            qchem_out_path  - The path to the output file to parse.

        Returns:
            The energy result in the output file.
        """

        # find the energy inside the qchem file output
        with open(qchem_out_path) as qchem_out_file:
            for line in qchem_out_file:
                if line.find("Total energy in the final basis set = ") != -1:
                    return float(line[line.find("Total energy in the final basis set = ") + 39:])

        # if no line with "Total energy in the final basis set = " is found, raise an exception
        raise LibraryCallError("qchem", "energy calculation", "output file is of incorrect format", log_path=qchem_out_path)

    def optimize_geometry(self, molecule, model):
        """
        Optimizes the given input geometry with a provided model.

        Args:
            molecule        - The Molecule to perform the optimization on.
            model           - The Model to use for the optimization.

        Returns:
            (new optimized molecule, nergy of the optimized geometry, path to the log file)
        """

        # Check if the user has qchem isntalled

        print("Beginning geometry optimization using qchem of {} with {}/{}.".format(molecule.get_name(), model.get_method(), model.get_basis()))

        qchem_in_path = files.get_optimization_log_path(self.settings.get("files", "log_path"), molecule, model.get_method(), model.get_basis(), "in")

        self.create_input_file(qchem_in_path, molecule, model, "opt")

        qchem_out_path = files.get_optimization_log_path(self.settings.get("files", "log_path"), molecule, model.get_method(), model.get_basis(), "out")

        num_threads = self.settings.getint("qchem", "num_threads", 1)

        # make the qchem system call
        try:
            system.call("qchem", "-nt", str(num_threads), qchem_in_path, qchem_out_path)
        except CommandExecutionError as e:
            raise LibraryCallError("qchem", "optimize", str(e), log_path=qchem_out_path)

        geometry, energy = self.find_optimized_geometry_and_energy_in_optimization_output_file(qchem_out_path, molecule.get_num_atoms(), molecule.get_charge(), molecule.get_spin_multiplicity())

        if self.logging:
            print("Completed geometry optimization.")

        return geometry, energy, qchem_out_path
        
    def find_optimized_geometry_and_energy_in_optimization_output_file(self, qchem_out_path, num_atoms, charge, spin_multiplicity):
        """
        Parses the output file of a Qchem geometry optimization and retrieves the optimized
        geometry and optimized energy that resulted from the calculation.

        Args:
            qchem_out_path  - The path to the output file to parse.

        Returns:
            (optimized geometry Molecule, optimized energy) from the output file.
        """

        # parse the molecule's geometry out of the output file
        found = False
        # found is set to true when the keyword 'Final energy is' is found. If found is true, it then looks for the keyword
        # 'ATOM' to look for the optimized geometry
        qchem_out_string = ""
        with open(qchem_out_path, "r") as qchem_out_file:
            for line in qchem_out_file:
                if found:
                    if "ATOM" in line:
                        for atom_index in range(num_atoms):
                            qchem_out_string += " ".join(qchem_out_file.readline().split()[1:]) + "\n"
                        return Molecule().read_psi4_string("{} {}\n{}".format(charge,
                                spin_multiplicity, qchem_out_string)), energy
                elif "Final energy is" in line:
                    energy = float(line.split()[3])
                    found = True

        # if we didn't find a geometry to parse, raise an exception.
        raise LibraryCallError("qchem", "optimze", "output file is of incorrect format", log_path=qchem_out_path)

    def find_frequencies(self, molecule, model):
        """
        Performs a frequency calculation to find the normal modes, frequencies, and reduced masses of the molecule

        Args:
            molecule        - The Molecule to perform the frequency calculation on.
            model           - The Model to use for the claculation.

        Returns:
            (normal modes, frequencies, reduced masses, path to log file)
        """

        print("Beginning normal modes calculation of {} with {}/{}.".format(molecule.get_name(), model.get_method(),
            model.get_basis()))

        qchem_in_path = files.get_frequencies_log_path(self.settings.get("files", "log_path"), molecule, model.get_method(), model.get_basis(), "in")

        self.create_input_file(qchem_in_path, molecule, model, "freq")

        qchem_out_path = files.get_frequencies_log_path(self.settings.get("files", "log_path"), molecule, model.get_method(), model.get_basis(), "out")

        num_threads = self.settings.getint("qchem", "num_threads", 1)

        try:
            system.call("qchem", "-nt", str(num_threads), qchem_in_path, qchem_out_path)
        except CommandExecutionError as e:
            raise LibraryCallError("qchem", "frequency", str(e), log_path=qchem_out_path)


        normal_modes, frequencies, red_masses = self.find_normal_modes_frequencies_and_reduced_masses_in_frequency_output_file(qchem_out_path, molecule.get_num_atoms())

        print("Normal mode/frequency analysis complete. {} normal modes found".format(len(normal_modes)))

        self.check_neg_freqs(frequencies)

        return normal_modes, frequencies, red_masses, qchem_out_path

    def find_normal_modes_frequencies_and_reduced_masses_in_frequency_output_file(self, qchem_out_path, num_atoms):
        """
        Parses the output file of a Qchem frequency calculation and retrieves the normal modes, frequencies,
        and reduced masses

        Args:
            qchem_out_path  - The path to the output file to parse.

        Returns:
            (normal modes, frequencies, reduced masses) from the output file.
        """

        # now parse the qchem output file to read the normal modes, frequencies, and reduced masses
        with open(qchem_out_path, "r") as qchem_out_file:
            frequencies = []
            red_masses = []
            normal_modes = []
            while(True):
                line = qchem_out_file.readline()

                if line == "":
                    if frequencies == []:
                        raise LibraryCallError("qchem", "frequency", "process returned file of incorrect format")
                    break

                if "Mode:" in line:
                    num_modes = len(line.split()) - 1

                    # read frequencies from file and cast each to a float
                    frequencies += [float(freq) for freq in qchem_out_file.readline().split()[1:]]

                    # skip the Force Constant line
                    qchem_out_file.readline()

                    # read red masses from file and cast each to float
                    red_masses += [float(freq) for freq in qchem_out_file.readline().split()[2:]]

                    # skip IR Active, IR Intens, Raman Active and labels line
                    for i in range(4):
                        qchem_out_file.readline()

                    # initialize 3d list for normal modes
                    modes = [[[0, 0, 0] for i in range(num_atoms)] for k in range(num_modes)]

                    for atom in range(num_atoms):

                        # read the coordinates for the next atom for the next num_modes modes
                        modes_of_atom = [float(i) for i in qchem_out_file.readline().split()[1:]]

                        # loop over each mode
                        for index in range(num_modes):

                            # store coordinates in corresponding entries of modes list
                            modes[index][atom][0] = modes_of_atom[index * 3]
                            modes[index][atom][1] = modes_of_atom[index * 3 + 1]
                            modes[index][atom][2] = modes_of_atom[index * 3 + 2]

                    normal_modes += modes

        return normal_modes, frequencies, red_masses
