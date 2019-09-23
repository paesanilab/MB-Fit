# External package imports
import math
from random import Random

# Absolute module impots
from potential_fitting.utils import Quaternion
from potential_fitting.molecule import Molecule
from potential_fitting.utils import files, system
from .configuration_generator import ConfigurationGenerator
from potential_fitting.utils.distribution_function import LinearDistributionFunction, LogarithmicDistributionFunction, RandomDistributionFunction

class DistanceSamplingConfigurationGenerator(ConfigurationGenerator):
    """
    Implementation of ConfigurationGenerator that generates configurations using
    a random sampling over a variety of distances.
    """

    def __init__(self, settings_path, min_distance=1, max_distance=5, min_inter_distance=0.8, progression=False,
                 use_grid=False, step_size=0.5, num_attempts=100, logarithmic=False, distribution=None):
        """
        Constructs a new DistanceSamplingConfigurationGenerator.

        Args:
            settings_path       - Local path to '.ini' settings file with all relevant settings.
            min_distance        - The minimum distance between the centers of mass of the two monomers.
            max_distance        - The maximum distance between the centers of mass of the two monomers.
            min_inter_distance  - Minimum intermolecular distance is this times the sum of the van der Waals radii of two
                    atoms.
            progression         - If True, a smooth progression will be used over the chosen distribution. Otherwise,
                    random points on the distribution will be sampled.
            use_grid            - If True, then distance between the center of mass of the monomers will be on a grid.
                    Otherwise, it will be smooth.
            step_size           - Only used if use_grid is True, this is the step size of the grid in angstroms.
            num_attempts        - The number of attempts to generate a configuration at any given distance before giving
                    up and moving to the next distance.
            logarithmic         - If True, then a logarithmic progression is used to generate the configurations.
                    This means more configs are generated at lower distances.
            distribution        - An implementation of DistributionFunction. If specified, the logarithmic argument
                    is ignored and this distribution is used to choose the distances between configurations. Should
                    be implemented over the domain [0,1]. So the first config will have distance =
                    distribution.get_value(0) and the last config will have distance = distribution.get_value(1).

        Returns:
            A new DistanceSamplingConfigurationGenerator.
        """

        super(DistanceSamplingConfigurationGenerator, self).__init__(settings_path)

        self.min_inter_distance = min_inter_distance
        self.use_grid = use_grid
        if max_distance == min_distance:
            self.step_size = 1
        else:
            self.step_size = step_size / (max_distance - min_distance)
        self.num_attempts = num_attempts
        self.random = Random()

        if distribution is not None:
            self.distance_distribution=distribution
        elif logarithmic:
            self.distance_distribution = LogarithmicDistributionFunction(min_distance, max_distance, 0, 1)
        else:
            self.distance_distribution = LinearDistributionFunction.get_function_from_2_points(0, min_distance, 1, max_distance)

        if not progression:
            self.distance_distribution = RandomDistributionFunction(self.distance_distribution, self.random, 0, 1)

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

        # put both geometries in standard orientation if they are not already
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

        # if use_grid is false, set the step size to even space the configurations
        if not self.use_grid:
            step_size = 1 / num_configs
        else:
            step_size = self.step_size

        # how many steps the grid will have
        num_steps = math.floor(1 / step_size)

        # keeps track of how many total configurations have been generated
        total_configs = 0

        # construct a psuedo-random number generator
        self.random.seed(seed)

        system.format_print(
            "Beginning 2B configurations generation with smooth distribution. Will generate {} configs.".format(
                num_configs),
            bold=True, color=system.Color.YELLOW)


        # try up to 3 times to generate configs over the distribution.
        for cycle_index in range(0, 3):

            # loop over each step on our grid
            for step in range(num_steps):

                # loop over how many configs we want to generate at this step in the grid, which is equal
                #   to the number of configs remaining to be generated divided by the number of steps left.
                #   this ensures that unless a config at the last step is impossible, we will always have
                #   exactly num_configs configs.

                for config in range(math.ceil((num_configs - total_configs) / (num_steps - step))):

                    distance = self.distance_distribution.get_value(step * step_size)

                    # first select a random geometry for each monomer
                    molecule1 = self.random.choice(molecules1)
                    molecule2 = self.random.choice(molecules2)

                    try:
                        self.move_to_config(self.random, molecule1, molecule2, distance)

                    except RanOutOfAttemptsException:
                        # if we didn't find a valid configuration, skip this config
                        continue

                    mol = Molecule.read_xyz_direct(str(molecule1.get_num_atoms() + molecule2.get_num_atoms()) + "\n\n" + molecule1.to_xyz() + "\n" + molecule2.to_xyz())

                    yield mol

                    total_configs += 1

                    if total_configs % 100 == 0:
                        system.format_print("{} configs done...".format(total_configs),
                                            italics=True)

                    # if we have hit our target number of configs, return
                    if total_configs == num_configs:
                        system.format_print("Done! Generated {} configurations.".format(num_configs), bold=True,
                                            color=system.Color.GREEN)
                        return

        system.format_print("Done! Generated {} configurations.".format(total_configs), bold=True,
                            color=system.Color.GREEN)

        if total_configs < num_configs:
            system.format_print("Generated fewer than {} configs because it was too hard to generate configs over the given distribution.".format(num_configs), bold=True,
                                color=system.Color.GREEN)


class RanOutOfAttemptsException(Exception):
    """
    Used to check if move_to_config runs out of attempts
    """
    pass
