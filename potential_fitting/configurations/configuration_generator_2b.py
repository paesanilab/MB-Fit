# External package imports
import math
from random import Random

# Absolute module impots
from potential_fitting.utils import Quaternion
from potential_fitting.molecule import Molecule
from potential_fitting.utils import files, system
from .configuration_generator import ConfigurationGenerator

class DistanceSamplingConfigurationGenerator(ConfigurationGenerator):
    """
    Implementation of ConfigurationGenerator that generates configurations using
    a random sampling over a variety of distances.
    """

    def __init__(self, settings_path, min_distance=1, max_distance=5, min_inter_distance=0.8, progression=False,
                 use_grid=False, step_size=0.5, num_attempts=100, logarithmic=False):
        """
        Constructs a new DistanceSamplingConfigurationGenerator.

        Args:
            settings_path       - Local path to '.ini' settings file with all relevant settings.
            min_distance        - The minimum distance between the centers of mass of the two monomers.
            max_distance        - The maximum distance between the centers of mass of the two monomers.
            min_inter_distance  - Minimum intermolecular distance is this times the sum of the van der walls radii of two
                    atoms.
            progression         - If True, a smooth progression will be used from min_distance to max_distance.
                    Otherwise, a random progression will be used.
            use_grid            - If True, then distance between the center of mass of the monomers will be on a grid.
                    Otherwise, it will be smooth.
            step_size           - Only used if use_grid is True, this is the step size of the grid in angstroms.
            num_attempts        - The number of attempts to generate a configuration at any given distance before giving
                    up and moving to the next distance.
            logarithmic         - If True, then a logarithmic progression is used to generate the configurations.
                    This means more configs are generated at lower distances.

        Returns:
            A new DistanceSamplingConfigurationGenerator.
        """

        super(DistanceSamplingConfigurationGenerator, self).__init__(settings_path)

        self.min_distance = min_distance
        self.max_distance = max_distance
        self.min_inter_distance = min_inter_distance
        self.progression = progression
        self.use_grid = use_grid
        self.step_size = step_size
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

    def move_to_config(self, random, molecule1, molecule2, distance):
        """
        Moves the molecules to a configuration with the given distance between their centers of
        mass. Raises RanOutOfAttemptsException if a configuration is failed to be found after a certain number of attempts.

        Args:
            random              - The Random object used to generate the configuration.
            molecule1           - Molecule object of the 1st monomer.
            molecule2           - Molecule object of the 2nd monomer.
            distance            - Distance between the centers of mass of the two molecules.

        Returns:
            None.
        """

        # put both geoemtries in standard orientation if they are not already
        molecule1.move_to_center_of_mass()
        molecule2.move_to_center_of_mass()
        molecule1.rotate_on_principal_axes()
        molecule2.rotate_on_principal_axes()

        # move the 2nd molecule away from the first
        molecule2.translate(distance, 0, 0)

        for attempt in range(self.num_attempts):

            # rotate each molecule a random amount
            molecule1.rotate(Quaternion.get_random_rotation_quaternion(random), 0, 0, 0)
            molecule2.rotate(Quaternion.get_random_rotation_quaternion(random), distance, 0, 0)

            # setting a flag variable to keep track of a valid configuration
            flag = True

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
                    and molecule_lists[1] is the same for the second molecule.
            num_configs     - The number of configurations to generate.
            seed            - Seed for the random number generator. The same seed will yield the same configurations
                    when all else is held equal.

        Yields:
            Molecule objects containing the new configurations.
        """

        if seed is None:
            seed = self.get_rand_seed()

        molecules1 = molecule_lists[0]
        molecules2 = molecule_lists[1]

        if self.progression:
            yield from self.generate_configurations_smooth(molecules1, molecules2, num_configs, seed)
        else:
            yield from self.generate_configurations_random(molecules1, molecules2, num_configs, seed)

    def generate_configurations_smooth(self, molecules1, molecules2, num_configs, seed):
        """
        Generate a set of 2 body configurations of the two geometries based on a smooth
        progression.

        Args:
            molecules1      - List of Molecule object configurations for the first molecule that will be sampled
                    from to generate configurations.
            molecules2      - List of Molecule object configurations for the second molecule that will be sampled
                    from to generate configurations.
            num_configs     - The number of configurations to generate.
            seed            - Seed for the random number generator. The same seed will yield the same configurations
                    when all else is held equal.

        Yields:
            Molecule objects containing the new configurations.
        """

        # if use_grid is false, set the step size to even space the configurations
        if not self.use_grid:
            step_size = (self.max_distance - self.min_distance) / num_configs

        # how many steps the grid will have
        num_steps = math.floor((self.max_distance - self.min_distance) / step_size)

        # keeps track of how many total configurations have been generated
        total_configs = 0

        # construct a psuedo-random number generator
        random = Random(seed)

        system.format_print(
            "Beginning 2B configurations generation with smooth distribution. Will generate {} configs.".format(
                num_configs),
            bold=True, color=system.Color.YELLOW)

        # loop over each step on our grid
        for step in range(num_steps):
            # loop over how many configs we want to generate at this step in the grid, which is equal
            #   to the number of configs remaining to be generated divided by the number of steps left.
            #   this ensures that unless a config at the last step is impossible, we will always have
            #   exactly num_configs configs.
            for config in range(math.ceil((num_configs - total_configs) / (num_steps - step))):

                # first select a random geometry for each monomer
                molecule1 = random.choice(molecules1)
                molecule2 = random.choice(molecules2)

                try:
                    # move the molecules to a valid configuration, making 5 attempts
                    if self.logarithmic:
                        self.move_to_config(random, molecule1, molecule2,
                                            self.logarithmic_progression(self.min_distance, self.max_distance, num_configs,
                                                               num_configs * step / (num_steps)))
                    else:
                        self.move_to_config(random, molecule1, molecule2, self.min_distance + step * step_size)

                except RanOutOfAttemptsException:
                    # if we didn't find a valid configuration, skip this config
                    continue

                mol = Molecule.read_xyz_direct(str(molecule1.get_num_atoms() + molecule2.get_num_atoms()) + "\n\n" + molecule1.to_xyz() + "\n" + molecule2.to_xyz())

                yield mol

                total_configs += 1

                if total_configs % 100 == 0:
                    system.format_print("{} configs done...".format(num_configs - total_configs),
                                        italics=True)

                # if we have hit our target number of configs, return
                if total_configs == num_configs:
                    system.format_print("Done! Generated {} configurations.".format(num_configs), bold=True,
                                        color=system.Color.GREEN)
                    return

    def generate_configurations_random(self, molecules1, molecules2, num_configs, seed):
        """
        Generate a set of 2 body configurations of the two geometries based on a random
        progression.

        Args:
            molecules1      - List of Molecule object configurations for the first molecule that will be sampled
                    from to generate configurations.
            molecules2      - List of Molecule object configurations for the second molecule that will be sampled
                    from to generate configurations.
            num_configs     - The number of configurations to generate.
            seed            - Seed for the random number generator. The same seed will yield the same configurations
                    when all else is held equal.

        Yields:
            Molecule objects containing the new configurations.
        """

        # construct a psuedo-random number generator
        random = Random(seed)

        # setting the total number of configs
        total_configs = num_configs

        system.format_print(
            "Beginning 2B configurations generation with random distribution. Will generate {} configs.".format(
                num_configs),
            bold=True, color=system.Color.YELLOW)

        while total_configs > 0:

            # random distance between min_distance and max_distance

            if self.logarithmic:
                random_distance = self.logarithmic_progression(self.min_distance, self.max_distance, num_configs,
                                                          random.uniform(0, num_configs - 1))
            else:
                random_distance = random.uniform(self.min_distance, self.max_distance)

            # getting the molecules
            molecule1 = random.choice(molecules1)
            molecule2 = random.choice(molecules2)

            # generating one confiugration at that random distance

            try:
                self.move_to_config(random, molecule1, molecule2, random_distance)
            except RanOutOfAttemptsException:
                # if we didn't find a valid configuration, skip this config
                continue

            mol = Molecule.read_xyz_direct(str(molecule1.get_num_atoms() + molecule2.get_num_atoms()) + "\n\n" + molecule1.to_xyz() + "\n" + molecule2.to_xyz())

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
