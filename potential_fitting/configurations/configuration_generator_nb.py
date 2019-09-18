# External package imports
import math
from random import Random, randint

# Absolute module impots
from potential_fitting.utils import Quaternion
from potential_fitting.molecule import Molecule, xyz_to_molecules
from potential_fitting.utils import files, SettingsReader, system

from .configuration_generator import ConfigurationGenerator

class RandomSamplingConfigurationGenerator(ConfigurationGenerator):
    """
    Implementation of ConfigurationGenerator that generates configurations by randomly placing molecules
    within a sphere of a given radius.
    """

    def __init__(self, settings_path, radius=10, min_inter_distance=0.8, num_attempts=100, logarithmic=False):
        """
        Constructs a new RandomSamplingConfigurationGenerator.

        Args:
            settings_path       - Local path to '.ini' settings file with all relevant settings.
            radius              - Radius of the sphere to place molecules in.
            min_inter_distance  - Minimum intermolecular distance is this times the sum of the van der walls radii of two
                    atoms.
            num_attempts        - The number of attempts to generate a configuration at any given distance before giving
                    up and moving to the next distance.
            logarithmic         - If True, then a logarithmic progression is used to generate the configurations.
                    This means more configs are generated closer to the center of the sphere.

        Returns:
            A new RandomSamplingConfigurationGenerator.
        """

        super(ConfigurationGenerator, self).__init__(settings_path)

        self.radius = radius
        self.min_inter_distance = min_inter_distance
        self.num_attempts = num_attempts
        self.logarithmic = logarithmic

    def logarithmic_progression(self, min, max, num_configs, config_num):
        """
        Gives a distance of a point on a logarithmic distribution.

        Args:
            min                 - The minimun distance.
            max                 - The maximum distance.
            num_configs         - The total number of configs in this progression.
            config_num          - The number of the current config, in range [0, num_configs - 1]
                    Config 0 will be at min distance while config num_configs - 1 will be at max distance.
        """

        dx = (math.log(max) - math.log(min)) / (num_configs - 1)
        return math.e ** (math.log(min) + config_num * dx)

    def move_to_config(self, random, molecules, distances):
        """
        Moves the molecules to a configuration with the given distance between their centers of
        mass and the point 0,0,0. Raises RanOutOfAttemptsException if a configuration is failed to be found after a certain number of attempts.

        Args:
            random              - The Random object used to generate the configuration.
            molecules           - The molecules to place in the configuration.
            distances           - Distance of each monomer from 0,0,0.

        Returns:
            None.
        """

        for attempt in range(self.attempts):

            for molecule, distance in zip(molecules, distances):
                molecule.move_to_center_of_mass()
                molecule.rotate_on_principal_axes()

                molecule.rotate(Quaternion.get_random_rotation_quaternion(random), 0, 0, 0)
                molecule.translate(*Quaternion.get_random_rotation_quaternion(random).rotate(distance, 0, 0))

            # setting a flag variable to keep track of a valid configuration
            flag = True

            for index1, molecule1 in enumerate(molecules):
                for index2, molecule2 in enumerate(molecules[index1 + 1:]):
                    # loop through atoms in molecule 1
                    for atom1 in molecule1.get_atoms():
                        # loop through atoms in molecule 2
                        for atom2 in molecule2.get_atoms():
                            # get the sum of the venderwall radii of the pair
                            sum_vdw_distance = atom1.get_vdw_radius() + atom2.get_vdw_radius()
                            # comparing with our factor (default = 0.8)
                            if atom1.distance(atom2) < self.min_inter_distance * sum_vdw_distance:
                                flag = False

            # Returning only valid configurations.
            if flag:
                return

        # if we run out of attempts without generating a valid configuration, raise an exception
        raise RanOutOfAttemptsException

    def generate_configurations(self, molecule_lists, num_configs, seed=None):
        """
        Generates configurations of the given molecule.

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

        if seed is None:
            seed = self.get_rand_seed()

        # construct a psuedo-random number generator
        random = Random(seed)

        # setting the total number of configs
        total_configs = num_configs

        system.format_print(
            "Beginning nB configurations generation. Will generate {} configs.".format(num_configs), bold=True,
            color=system.Color.GREEN)

        while total_configs > 0:

            # random distance between min_distance and max_distance

            if self.logarithmic:
                distances = [self.logarithmic_progression(0, self.radius, num_configs,
                                                     random.uniform(0, num_configs - 1)) for molecules in
                             molecule_lists]
            else:
                distances = [random.uniform(0, self.radius) for molecules in molecule_lists]

            molecules = [random.choice(molecules_list) for molecules_list in molecule_lists]

            # generating one confiugration at that random distance

            try:
                self.move_to_config(random, molecules, distances, self.min_inter_distance, self.num_attempts)
            except RanOutOfAttemptsException:
                # if we didn't find a valid configuration, skip this config
                continue

            mol = Molecule([molecule.get_fragments()[0] for molecule in molecules])
            yield mol

            # decrementing required number of configs
            total_configs -= 1

            if (num_configs - total_configs) % 100 == 0:
                system.format_print("{} configs done...".format(num_configs - total_configs), italics=True)

        # if we have hit our target number of configs, return
        system.format_print("Done! Generated {} configurations.".format(num_configs), bold=True,
                            color=system.Color.GREEN)
        return


class RanOutOfAttemptsException(Exception):
    """
    Used to check if move_to_config runs out of attempts
    """
    pass
