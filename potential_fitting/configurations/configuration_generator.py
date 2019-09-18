# external package imports
from random import randint

# absolute module imports
from potential_fitting.utils import SettingsReader, system

class ConfigurationGenerator(object):
    """
    Abstract class that all configuration generators will extend.
    """

    def __init__(self, settings_path):
        """
        Constructs a new ConfigurationGenerator.

        Args:
            settings_path       - Local path to '.ini' settings file with all relevant settings.

        Returns:
            A new ConfigurationGenerator.
        """

        self.settings = SettingsReader(settings_path)

    def generate_configurations(self, molecule_lists, num_configs, seed=None):
        """
        Generates Configurations of one or more molecules.

        Args:
            molecule_lists  - List of lists of molecules to generate configurations from such that molecule_lists[0]
                    is a list of all configurations to use in the generation of configurations for the first molecule
                    and so on.
            num_configs     - The number of configurations to generate.
            seed            - Seed for the random number generator. The same seed will yield the same configurations
                    when all else is held equal.

        Yields:
            Molecule objects containing the new configurations.
        """

        raise NotImplementedError

    def get_rand_seed(self):
        """
        Gets a new random seed and logs a message to the console.

        Args:
            None.

        Returns:
            A new random seed in the range [-100000, 100000].
        """

        seed = randint(-100000, 100000)
        system.format_print("No seed given, randomly using seed {}.".format(seed),
                            italics=True)
        return seed
