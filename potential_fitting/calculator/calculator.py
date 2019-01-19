# absolute module imports
from potential_fitting.utils import SettingsReader
from potential_fitting.molecule import Molecule

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

        raise NotImplementedError

    def optimize_geometry(self, molecule, model):
        """
        Optimizes the given input geometry with a provided model.

        Args:
            molecule        - The Molecule to perform the optimization on.
            model           - The Model to use for the optimization.

        Returns:
            (new optimized molecule, nergy of the optimized geometry, path to the log file)
        """

        raise NotImplementedError

    def find_frequencies(self, molecule, model):
        """
        Performs a frequency calculation to find the normal modes, frequencies, and reduced masses of the molecule

        Args:
            molecule        - The Molecule to perform the frequency calculation on.
            model           - The Model to use for the claculation.

        Returns:
            (normal modes, frequencies, reduced masses, path to log file)
        """

        raise NotImplementedError

    def check_neg_freqs(self, frequencies):
        number_of_negative_frequencies = 0
        for frequency in frequencies:
            if frequency < 0:
                number_of_negative_frequencues += 1

        if number_of_negative_frequencies = 1:
            print("Single negative frequency detected. This means the inputted geometry is probably a transition state.")
        elif number_of_negative_frequencies > 1:
            print("Multiple ({}) negative frequencies detected. Proceed with caution.")
