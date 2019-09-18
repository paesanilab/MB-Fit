# external package imports
from random import randint

# absolute module imports
from potential_fitting.utils import SettingsReader, system, files
from potential_fitting.molecule import xyz_to_molecules

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

    @staticmethod
    def generate_configs_from_file_to_file(geo_paths, out_path, config_generator, num_configs, seed=None):
        """
        This function reads geometries from filepaths, feeds them into a config_generator, and then writes
        the output to another file.

        Args:
            geo_paths       - List of local paths to '.xyz' files containing geometries to generate configurations with.
            out_path        - Local path to '.xyz' file to write configurations to.
            config_generator - Implementation of ConfigurationGenerator to use to generate the configurations.
            num_configs     - The number of configurations to generate.
            seed            - Seed given to config_generator. The same seed will produce the same configurations when
                    all else is held equal.
        """

        molecule_lists = [xyz_to_molecules(path) for path in geo_paths]

        molecules = config_generator.generate_configurations(molecule_lists, num_configs, seed=seed)

        out_path = files.init_file(out_path)

        with open(out_path, "w") as out_file:
            for molecule_index, molecule in enumerate(molecules):
                out_file.write("{}\n{}\n{}\n".format(molecule.get_num_atoms(), molecule_index, molecule.to_xyz()))
