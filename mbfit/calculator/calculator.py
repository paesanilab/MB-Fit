# absolute module imports
from mbfit.utils import SettingsReader, files
from mbfit.exceptions import NoSuchLibraryError
from mbfit.molecule import Molecule

def generate_input(molecule, fragment_indicies, model, cp, settings):
    """
    Generates input for energy of a subset of the fragments of a molecule.
    Args:
        molecule            - The Molecule object to calculate the energy of.
        fragment_indicies   - List of indicies of fragments to include in the calculation.
        model               - The model to use for this calculation, should be specified as method/basis.
        cp                  - Whether to use counterpoise correction, True means that fragments not included in
                fragment_indicies will be added as ghost atoms to the calculation.
        settings            - SettingsReader object with all relevent settings.

    Returns:
        Nothing
    """
    
    code = settings.get("energy_calculator", "code")
    method, basis = model.split("/")
    if code == "qchem":
        suffix = ""
        for frag in fragment_indicies:
            suffix +=  "_" + str(frag)
        suffix += ".in"
        # file to write qchem input in
        qchem_in_path = files.get_energy_log_path(settings.get("files", "log_path"), molecule, method, basis, cp, suffix)
        qchem_input = files.get_qchem_input_string(molecule, fragment_indicies, model, cp, settings)
        with open(qchem_in_path, "w") as qchem_in_file:
            qchem_in_file.write(qchem_input)

def run_calculation(molecule, fragment_indicies, model, cp, settings):
    """
    Calculates the energy of a subset of the fragments of a molecule. Assumes generate_input has been run.

    Args:
        molecule            - The Molecule object to calculate the energy of.
        fragment_indicies   - List of indicies of fragments to include in the calculation.
        model               - The model to use for this calculation, should be specified as method/basis.
        cp                  - Whether to use counterpoise correction, True means that fragments not included in
                fragment_indicies will be added as ghost atoms to the calculation.
        settings            - SettingsReader object with all relevent settings.

    Returns:
        Nothing
    """
    code = settings.get("energy_calculator", "code")
    method, basis = model.split("/")
    if code == "qchem":
        suffix = ""
        for frag in fragment_indicies:
            suffix +=  "_" + str(frag)
        
        qchem_in_path = files.get_energy_log_path(settings.get("files", "log_path"), molecule, method, basis, cp, suffix + ".in")
        qchem_out_path = files.get_energy_log_path(settings.get("files", "log_path"), molecule, method, basis, cp, suffix + ".out")

        # get number of threads
        num_threads = settings.getint("qchem", "num_threads")
        # perform system call to run qchem
        try:
            system.call("qchem", "-nt", str(num_threads), qchem_in_path, qchem_out_path)
        except CommandExecutionError as e:
            raise LibraryCallError("qchem", "energy calculation", "process returned non-zero exit code") from e

def retrieve_energy(molecule, fragment_indicies, model, cp, settings):
    """
    Retrieves the energy of a subset of the fragments of a molecule. Assumes run_calculation has been run.

    Args:
        molecule            - The Molecule object to calculate the energy of.
        fragment_indicies   - List of indicies of fragments to include in the calculation.
        model               - The model to use for this calculation, should be specified as method/basis.
        cp                  - Whether to use counterpoise correction, True means that fragments not included in
                fragment_indicies will be added as ghost atoms to the calculation.
        settings            - SettingsReader object with all relevent settings.

    Returns:
        The calculated energy of the subset of fragments in this molecule.
    """
    code = settings.get("energy_calculator", "code")
    method, basis = model.split("/")
    if code == "qchem":
        suffix = ""
        for frag in fragment_indicies:
            suffix += "_" + str(frag)

        qchem_in_path = files.get_energy_log_path(settings.get("files", "log_path"), molecule, method, basis, cp, suffix + ".in")
        qchem_out_path = files.get_energy_log_path(settings.get("files", "log_path"), molecule, method, basis, cp, suffix + ".out")

        # find the energy inside the qchem file output
        with open(qchem_out_path) as qchem_out_file:
            for line in qchem_out_file:
                if line.find("Total energy in the final basis set = ") != -1:
                    return float(line[line.find("Total energy in the final basis set = ") + 39:])
    
        # if no line with "Total energy in the final basis set = " is found, raise an exception
        raise LibraryCallError("qchem", "energy calculation", "process returned file of incorrect format")

        
class Calculator:
    def __init__(self, settings_path, logging = True):
        """
        Constructor for a Calculator.

        Args:
            settings_path   - Path to the settings.ini
            logging         - Whether this calculator should output logging messages.
                Default: True

        Returns:
            None
        """

        self.settings = SettingsReader(settings_path)
        self.logging = logging

    def setLogging(self, logging):
        self.logging = logging

    def is_installed(self):
        raise NotImplementedError

    def is_valid_model(self, model):
        raise NotImplementedError
    
    def calculate_energy(self, molecule, model, fragment_indicies, qm_options={}):
        """
        Calculates the energy of a subset of the fragments of a molecule with a provided model.

        Fragments not included in fragment_indicies will be included as ghost atoms if model's
        cp is enabled.

        Args:
            molecule        - The Molecule to perform the calculation on.
            fragment_indicies - List of the indicies of the fragments to include in the calculation.
            model           - The Model to use for the calculation.
            qm_options       - Dictionary of extra arguments to be passed to the QM code doing the calculation.

        Returns:
            (calculated energy, path to the log file)
        """

        raise NotImplementedError

    def optimize_geometry(self, molecule, model, qm_options={}):
        """
        Optimizes the given input geometry with a provided model.

        Args:
            molecule        - The Molecule to perform the optimization on.
            model           - The Model to use for the optimization.
            qm_options       - Dictionary of extra arguments to be passed to the QM code doing the calculation.

        Returns:
            (new optimized molecule, energy of the optimized geometry, path to the log file)
        """

        raise NotImplementedError

    def calculate_frequencies(self, molecule, model, qm_options={}):
        """
        Performs a frequency calculation to find the normal modes, frequencies, and reduced masses of the molecule

        Args:
            molecule        - The Molecule to perform the frequency calculation on.
            model           - The Model to use for the claculation.
            qm_options       - Dictionary of extra arguments to be passed to the QM code doing the calculation.

        Returns:
            (normal modes, frequencies, reduced masses, path to log file)
        """

        raise NotImplementedError

    def check_neg_freqs(self, frequencies):
        number_of_negative_frequencies = 0
        for frequency in frequencies:
            if frequency < 0:
                number_of_negative_frequencies += 1

        if number_of_negative_frequencies == 1:
            print("Single negative frequency detected. This means the inputted geometry is probably a transition state.")
        elif number_of_negative_frequencies > 1:
            print("Multiple ({}) negative frequencies detected. Proceed with caution.".format(number_of_negative_frequencies))

        return number_of_negative_frequencies
