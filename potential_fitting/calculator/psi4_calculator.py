# absolute module imports
from from potential_fitting.utils import files
from potential_fitting.molecule import Molecule
from potential_fitting.exceptions import (LibraryNotAvailableError, LibraryCallError)

try:
    import psi4
    from psi4.driver.qcdb.exceptions import QcdbException
except ImportError:
    has_psi4 = False
else:
    has_psi4 = True

class Psi4Calclator(Calculator):
    
    def __init__(self, settings_path, logging = True):
        """
        Constructor for a Psi4Calculator.

        A Psi4Calculator uses Psi4 for all of its calculations.

        Args:
            settings_path   - Path to the settings.ini
            logging         - Whether this calculator should output logging messages.
                Default: True

        Returns:
            None
        """

        super(Psi4Calculator, self).__init__(settings_path, logging)

        if not has_psi4:
            raise LibraryNotAvailableError("psi4")

    def initialize_calculation(self, log_path):
        """
        Sets the log file, number of threads, and memory for psi4 to use.

        Args:
            log_path        - The path to the log file to set.

        Returns:
            None.
        """

        # set the log file
        psi4.core.set_output_file(log_path, False)

        # set the number of threads to use to compute based on settings
        psi4.set_num_threads(self.settings.getint("psi4", "num_threads", 1))

        # set the amount of memory to use based on settings
        psi4.set_memory(self.settings.get("psi4", "memory", "1000 MB"))

    def get_psi4_molecule(self, molecule, model, fragment_indicies = None):
        """
        Converts a provided PotentialFitting Molecule into a Psi4 Molecule needed for Psi4 calculations.

        Args:
            molecule        - The PotentialFitting Molecule to convert.
            model           - The Model to check if cp is enabled or not. If model has cp enabled, then fragments
                    not included in fragment_indicies will be included as ghost atoms.
            fragment_indicies - If not None, then use only these fragments will be included when translating to a Psi4 Molecule.
                Default: None

        Returns:
            The Psi4 Molecule converted from the provided PotentialFitting Molecule
        """
        
        # Creats the psi4 input string of the molecule by combining the xyz file output with an additional line containing
        # charge and spin multiplicity
        psi4_string = "{}\n{} {}".format(molecule.to_xyz(fragment_indicies, model.get_cp()), molecule.get_charge(fragment_indicies),
                molecule.get_spin_multiplicity(fragment_indicies))

        try:
            # Constructs the psi4 Molecule from the string representation by making a call to a psi4 library function
            return psi4.geometry(psi4_string)
        except RuntimeError as e:
            raise LibraryCallError("psi4", "geometry", str(e))

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
        
        # file to write logs from psi4 calculation
        log_path = files.get_energy_log_path(self.settings.get("files", "log_path"), molecule, model.get_method(), model.get_basis(), model.get_cp(), "out")
        
        self.initialize_calculation(log_path)

        psi4_mol = self.get_psi4_molecule(molecule, model, fragment_indicies)
        
        try:
            # Perform library call to calculate energy of molecule
            if self.logging:
                print("Successfully completed energy calculation.")

            return psi4.energy(model.get_method() + "/" + model.get_basis(), molecule=psi4_mol), log_path

        except QcdbException as e:
            raise LibraryCallError("psi4", "energy", str(e))

    def optimize_geometry(self, molecule, model):
        """
        Optimizes the given input geometry with a provided model.

        Args:
            molecule        - The Molecule to perform the optimization on.
            model           - The Model to use for the optimization.

        Returns:
            (new optimized molecule, nergy of the optimized geometry, path to the log file)
        """

        if self.logging:
            print("Beginning geometry optimization of {} with model {}/{}".format(molecule.get_name(), model.get_method(), model.get_basis()))

        log_path = files.get_optimization_log_path(self.settings.get("files", "log_path"), molecule, model.get_method(), model.get_basis(), "log")

        self.initialize_calculation(log_path)

        psi4_mol = self.get_psi4_molecule(molecule, model)

        # use psi4 to optimize the geometry
        try:
            energy = psi4.optimize("{}/{}".format(model.get_method(), model.get_basis()), molecule=psi4_mol)
        except QcdbException as e:
            raise LibraryCallError("psi4", "optimize", str(e))

        if self.logging:
            print("Completed geometry optimization.")

        return Molecule().read_psi4_string(psi4_mol.save_string_xyz()), energy, log_path

    def find_frequencies(self, molecule, model):
        """
        Performs a frequency calculation to find the normal modes, frequencies, and reduced masses of the molecule

        Args:
            molecule        - The Molecule to perform the frequency calculation on.
            model           - The Model to use for the claculation.

        Returns:
            (normal modes, frequencies, reduced masses, path to log file)
        """

        print("Beginning normal modes calculation of {} with {}/{}.".format(molecule.get_name(), method, basis))

        log_path = files.get_frequencies_log_path(settings.get("files", "log_path"), molecule, method, basis, "log")

        self.initialize_calculation(log_path)

        psi4_mol = self.get_psi4_molecule(molecule, model)

        # use psi4 to perform a frequency calculation
        try:
            total_energy, wavefunction = psi4.frequency("{}/{}".format(method, basis), molecule=psi4_mol, return_wfn=True)
        except QcdbException as e:
            raise LibraryCallError("psi4", "frequency", str(e))

        # retrieve the normal modes, frequencies, and reduced masses from the psi4 output object
        vibration_info = psi4.qcdb.vib.filter_nonvib(wavefunction.frequency_analysis)

        normal_modes, frequencies, red_masses = self.parse_psi4_output(vibration_info)

        print("Normal mode/frequency analysis complete. {} normal modes found.".format(num_modes))

        self.check_neg_freqs(frequencies)

        return normal_modes, frequencies, red_masses, log_path

    def parse_frequencies_output_object(self, vibration_info):
        """
        Parses the Psi4 output of a frequency calculation into lists of
        normal_modes, frequencies, and reduced masses

        Args:
            vibration_info  - The Psi4 output of a frequency calculation.

        Returns:
            (normal modes, frequencies, reduced masses)
        """

        frequencies = wavefunction.frequencies().to_array()
        red_masses = vib_info["mu"].data
        
        normal_modes_raw = vib_info['x'].data
        normal_modes = []
        
        num_modes = len(frequencies)
        num_atoms = molecule.get_num_atoms()

        for normal_mode_index in range(num_modes):
            normal_mode = [[0, 0, 0] for atom_index in range(num_atoms)]
            
            for atom_index in range(num_atoms):
                normal_mode[atom_index][0] = float(normal_modes_raw[3 * atom_index + 0][normal_mode_index])
                normal_mode[atom_index][1] = float(normal_modes_raw[3 * atom_index + 1][normal_mode_index])
                normal_mode[atom_index][2] = float(normal_modes_raw[3 * atom_index + 2][normal_mode_index])
            
            normal_modes.append(normal_mode)
        
        return normal_modes, frequencies, red_masses