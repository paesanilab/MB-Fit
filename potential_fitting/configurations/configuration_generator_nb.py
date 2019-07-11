# External package imports
import math
from random import Random, randint

# Absolute module impots
from potential_fitting.utils import Quaternion
from potential_fitting.molecule import Molecule, xyz_to_molecules
from potential_fitting.utils import files, SettingsReader

def generate_configurations(settings_path, number_of_configs, config_path, *geo_paths, radius = 10, min_inter_distance=0.8, num_attempts=100,
                                      seed=None, logarithmic = False):
    """
    Generates a set of n body configurations by randomly placing monomer geometries in a sphere.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        number_of_configs   - Number of configurations to generate.
        config_path         - Local path to the file to write the configurations.
        geo_paths           - Paths of all the geometries to make the nb configurations from. Can be single optimized
                geometry or a set of distorted geometries. Will take random configs from these files to make the
                nb configs.
        radius              - Radius of the sphere monomers are placed within.
        min_inter_distance  - Minimum intermolecular distance is this times the sum of the van der walls radii of two
                atoms.
        num_attempts        - The number of attempts to generate a configuration at any given distance before giving
                up and moving to the next distance.
        seed                - Seed to use, the same seed will give the same configurations.
        logarithmic         - If True, will use logarithmic progression to chose distances from center of sphere
                for monomers.

    Returns:
        None
    """

    settings = SettingsReader(settings_path)

    if seed is None:
        seed = randint(-100000, 100000)

    # parse the molecules from the input xyz files
    molecules_lists = [xyz_to_molecules(geo_path) for geo_path in geo_paths]

    # construct a psuedo-random number generator
    random = Random(seed)

    # setting the total number of configs
    total_configs = number_of_configs

    print("Beginning nB configurations generation.")

    # open the config file to write to
    with open(files.init_file(config_path, files.OverwriteMethod.get_from_settings(settings)), "w") as config_file:

        while total_configs > 0:

            # random distance between min_distance and max_distance

            if logarithmic:
                distances = [logarithmic_progression(0, radius, number_of_configs,
                                                          random.uniform(0, number_of_configs - 1)) for molecules in molecules_lists]
            else:
                distances = [random.uniform(0, radius) for molecules in molecules_lists]

            molecules = [random.choice(molecules_list) for molecules_list in molecules_lists]

            # generating one confiugration at that random distance

            try:
                move_to_config(random, molecules, distances, min_inter_distance, num_attempts)
            except RanOutOfAttemptsException:
                # if we didn't find a valid configuration, skip this config
                continue
            # write total number of atoms to config file
            config_file.write("{}\n".format(sum([molecule.get_num_atoms() for molecule in molecules])))

            # in the comment line, write how many configs have been generated before this one
            config_file.write("{}\n".format(number_of_configs - total_configs))

            # write the xyz of each monomer to the config file
            for molecule in molecules:
                config_file.write("{}\n".format(molecule.to_xyz()))

            # decrementing required number of configs
            total_configs -= 1

            if (number_of_configs - total_configs) % 100 == 0:
                print("{} configs done...".format(number_of_configs - total_configs))

    # if we have hit our target number of configs, return
    print("Done! Generated {} configurations".format(number_of_configs))
    return


def logarithmic_progression(min, max, num_configs, config_num):
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


def move_to_config(random, molecules, distances, min_inter_distance, attempts):
    """
    Moves the molecules to a configuration with the given distance between their centers of
    mass and the point 0,0,0. Raises RanOutOfAttemptsException if a configuration is failed to be found after a certain number of attempts.

    Args:
        random              - The Random object used to generate the configuration.
        molecules           - The molecules to place in the configuration.
        distances           - Distance of each monomer from 0,0,0.
        min_inter_distance  - Minimum intermolecular distance is this times the sum of the atoms van der walls radii.
        attempts            - The number of configurations to try before giving up and raising the Exception.

    Returns:
        None.
    """

    for attempt in range(attempts):

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
                        if atom1.distance(atom2) < min_inter_distance * sum_vdw_distance:
                            flag = False

        # Returning only valid configurations.
        if flag:
            return

    # if we run out of attempts without generating a valid configuration, raise an exception
    raise RanOutOfAttemptsException


class RanOutOfAttemptsException(Exception):
    """
    Used to check if move_to_config runs out of attempts
    """
    pass
