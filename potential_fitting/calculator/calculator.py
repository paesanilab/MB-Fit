# external package imports
import os, subprocess

# absolute module imports
from potential_fitting.utils import SettingsReader, constants, files, system
from potential_fitting.exceptions import (LibraryNotAvailableError, LibraryCallError, NoSuchLibraryError,
        ConfigMissingSectionError, ConfigMissingPropertyError, CommandExecutionError)

has_psi4 = True

try:
    from psi4.driver.qcdb.exceptions import QcdbException
    import psi4
except ImportError:
    has_psi4 = False

class Calculator:

    def __init__(self, settings_path):
        """
        Constructor for a Calculator.

        Args:
            settings_path   - Path to the settings.ini

        Returns:
            None
        """

        self.settings = SettingsReader(settings_path)
    
    def calculate_energy(self, molecule, fragment_indicies, model):
        """
        Calculates the energy of a subset of the fragments of a molecule with a provided model.

        Fragments not included in fragment_indicies will be included as ghost atoms if model's
        cp is enables.

        Args:
            molecule        - The molecule to perform the calculation on.
            fragment_indicies - List of the indicies of the fragments to include in the calculation.
            model           - The model to use for the calculation.

        Returns:
            The calculated energy and the path to the log file
        """

        raise NotImplementedError

    def optimize_geometry(self, molecule, model):
        """
        Optimizes the given input geometry

        Args:
            molecule        - The molecule to perform the optimization on.
            model           - The model to use for the calculation.

        Returns:
            A new optimized molecule, the energy of the optimized geometry, and the path to the log file.
        """

        raise NotImplementedError

class Psi4Calclator(Calculator):
    
    def __init__(self, settings_path):
        """
        Constructor for a Psi4Calculator.

        Args:
            settings_path   - Path to the settings.ini

        Returns:
            None
        """

        super(Psi4Calculator, self).__init__(settings_path)

        if not has_psi4:
            raise LibraryNotAvailableError("psi4")

    def calculate_energy(self, molecule, fragment_indicies, model):
        """
        Calculates the energy of a subset of the fragments of a molecule with a provided model.

        Fragments not included in fragment_indicies will be included as ghost atoms if model's
        cp is enables.

        Args:
            molecule        - The molecule to perform the calculation on.
            fragment_indicies - List of the indicies of the fragments to include in the calculation.
            model           - The model to use for the calculation.

        Returns:
            The calculated energy and the path to the log file
        """
        
        # file to write logs from psi4 calculation
        log_path = files.get_energy_log_path(self.settings.get("files", "log_path"), molecule, model.get_method(), model.get_basis(), model.get_cp(), "out")
        
        # set the log file
        psi4.core.set_output_file(log_path, False)

        # Creats the psi4 input string of the molecule by combining the xyz file output with an additional line containing
        # charge and spin multiplicity
        psi4_string = "{}\n{} {}".format(molecule.to_xyz(fragment_indicies, model.get_cp()), molecule.get_charge(fragment_indicies),
                molecule.get_spin_multiplicity(fragment_indicies))

        try:
            # Constructs the psi4 Molecule from the string representation by making a call to a psi4 library function
            psi4_mol = psi4.geometry(psi4_string)
        except RuntimeError as e:
            raise LibraryCallError("psi4", "geometry", str(e))

        # Set the number of threads to use to compute based on settings
        psi4.set_num_threads(settings.getint("psi4", "num_threads"))

        # Set the amount of memory to use based on settings
        psi4.set_memory(settings.get("psi4", "memory"))
        
        try:
            # Perform library call to calculate energy of molecule
            return psi4.energy(model.get_method() + "/" + model.get_basis(), molecule=psi4_mol), log_path
        except QcdbException as e:
            raise LibraryCallError("psi4", "energy", str(e))

    def optimize_geometry(self, molecule, mode):
        """
        Optimizes the given input geometry

        Args:
            molecule        - The molecule to perform the optimization on.
            model           - The model to use for the calculation.

        Returns:
            A new optimized molecule, the energy of the optimized geometry, and the path to the log file.
        """

        # TODO: FINISH THIS AND OTHER METHODS


class QchemCalculator(Calculator):

    def __init__(self, settings_path):
        """
        Constructor for a QchemCalculator.

        Args:
            settings_path   - Path to the settings.ini

        Returns:
            None
        """
        
        super(QchemCalculator, self).__init__(settings_path)

        try:
            system.call("which", "qchem")
        except CommandExecutionError:
            raise LibraryNotAvailableError("qchem")

    def calculate_energy(self, molecule, fragment_indicies, model):
        """
        Calculates the energy of a subset of the fragments of a molecule with a provided model.

        Fragments not included in fragment_indicies will be included as ghost atoms if model's
        cp is enables.

        Args:
            molecule        - The molecule to perform the calculation on.
            fragment_indicies - List of the indicies of the fragments to include in the calculation.
            model           - The model to use for the calculation.

        Returns:
            The calculated energy and the path to the log file
        """


        # initialize qchem input string
        qchem_input = "";
        
        # molecule format
        qchem_input += "$molecule\n"

        # charge and spin multiplicity
        qchem_input += "{} {}\n".format(molecule.get_charge(fragment_indicies),
                molecule.get_spin_multiplicity(fragment_indicies))

        # atoms in the molecule
        # might need to add whitespace before each line?
        qchem_input += molecule.to_xyz(fragment_indicies, model.get_cp()) + "\n"

        qchem_input += "$end\n"

        # Q-chem settings
        qchem_input += "$rem\n"

        qchem_input += "jobtype " + "sp" + "\n"
        qchem_input += "method " + model.get_method() + "\n"
        qchem_input += "basis " + model.get_basis() + "\n"

        try:
            qchem_input += "ecp " + settings.get("qchem", "ecp") + "\n"
        except (ConfigMissingSectionError, ConfigMissingPropertyError):
            pass

        qchem_input += "$end"

        # file to write qchem input in
        qchem_in_path = files.get_energy_log_path(settings.get("files", "log_path"), molecule, molecule.get_method(), molecule.get_basis(), molecule.get_cp(), "in")
        
        with open(qchem_in_path, "w") as qchem_in_file:
            qchem_in_file.write(qchem_input)
        
        # file to write qchem output in
        qchem_out_path = files.get_energy_log_path(settings.get("files", "log_path"), molecule, molecule.get_method(), molecule.get_basis(), molecule.get_cp(), "out")

        # get number of threads
        num_threads = settings.getint("qchem", "num_threads")

        # perform system call to run qchem
        try:
            system.call("qchem", "-nt", str(num_threads), qchem_in_path, qchem_out_path)
        except CommandExecutionError as e:
            raise LibraryCallError("qchem", "energy calculation", "process returned non-zero exit code") from e

        # find the energy inside the qchem file output
        with open(qchem_out_path) as qchem_out_file:
            for line in qchem_out_file:
                if line.find("Total energy in the final basis set = ") != -1:
                    return float(line[line.find("Total energy in the final basis set = ") + 39:]), qchem_out_path

        # if no line with "Total energy in the final basis set = " is found, raise an exception
        raise LibraryCallError("qchem", "energy calculation", "process returned file of incorrect format")