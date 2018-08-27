import sys, os
import math
from random import Random, randint

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../qm_mb_energy_calculator/src/")

from molecule import Molecule
from molecule_parser import xyz_to_molecules
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import settings_reader

def generate_configurations(geo1, geo2, number_of_configs, config_path, min_distance = 1, max_distance = 5, min_inter_distance = 1.2, use_grid = False, step_size = 0.5, seed = randint(-100000, 1000000)):
    """
    Generates a set of 2 body configurations of the two optimized geometries and outputs them to an xyz file

    Args:
        geo1        - the first optimized (or series of unoptimized) geometry
        geo2        - the second optimized (or series of unoptimized) geometry
        number_of_configs - generate this many configurations
        config_path - the path to the file to write the configurations
        min_distance - the minimum distance between the centers of mass of the two monomers
        max_distance - the maximum distance between the centers of mass of the two monomers
        min_inter_distance - the minimum distance between any two atoms from oposite monomers
        use_grid    - if True, then distance between the center of mass of the monomers will be on a grid, otherwise, it will be smooth
        step_size - only used if use_grid is True, this is the step size of the grid in angstroms

    Returns:
        None
    """

    # if use_grid is false, set the step size to even space the configurations
    if not use_grid:
        step_size = (max_distance - min_distance) / number_of_configs

    # how many steps the grid will have
    num_steps = math.floor((max_distance - min_distance) / step_size)

    # how many configurations to generate per step in the grid (will be 1 if use_grid is False)
    configs_per_distance = math.ceil(number_of_configs / num_steps)

    # parse the molecules from the input xyz files
    molecules1 = xyz_to_molecules(geo1)
    molecules2 = xyz_to_molecules(geo2)

    # keeps track of how many total configurations have been generated
    total_configs = 0

    # construct a psuedo-random number generator
    random = Random(seed)

    # open the config file to write to
    with open(config_path, "w") as config_file:

        # loop over each step on our grid
        for step in range(num_steps):

            # loop over how many configs we want to generate at this step in the grid, which is equal
            #   to the number of configs remaining to be generated divided by the number of steps left.
            #   this ensures that unless a config at the last step is impossible, we will always have
            #   exactly number_of_configs configs.
            for config in range(math.ceil((number_of_configs - total_configs) / (num_steps - step))):

                # first select a random geometry for each monomer
                molecule1 = random.choice(molecules1)
                molecule2 = random.choice(molecules2)

                try:
                    # move the molecules to a valid configuration, making 5 attempts
                    move_to_config(random, molecule1, molecule2, min_distance + step * step_size, min_inter_distance, 5)

                except RanOutOfAttemptsException:
                    # if we didn't find a valid configuration, skip this config
                    continue

                # write total number of atoms to config file
                config_file.write("{}\n".format(molecule1.get_num_atoms() + molecule2.get_num_atoms()))

                # in the comment line, write how many configs have been generated before this one
                config_file.write("{}\n".format(total_configs))

                # write the xyz of each monomer to the config file
                config_file.write("{}\n{}\n".format(molecule1.to_xyz(), molecule2.to_xyz()))
            
                total_configs += 1

                # if we have hit our target number of configs, return
                if total_configs == number_of_configs:
                    print("Generated {} configurations".format(total_configs))
                    return

    # if we did not hit our target number of configs, notify the user
    print("Generated {} configurations".format(total_configs))

def move_to_config(random, molecule1, molecule2, distance, min_inter_distance, attempts):
    """
    Moves the given molecules to a configuration with the given distance between their centers of
    mass. Raises RanOutOfAttemptsException if a configuration is failed to be found.

    Args:
        random      - the Random object used to generate the configuration
        molecule1   - the Molecule object of the 1st monomer
        molecule2   - the Molecule object of the 2nd monomer
        distance    - distance between the centers of mass of the two molecules
        min_inter_distance - minimum distance for any inter-molecular atomic distance
        attempts    - the number of configurations to try before giving up

    Returns:
        None

    """

    # put both geoemtries in standard orientation if they are not already
    molecule1.move_to_center_of_mass()
    molecule2.move_to_center_of_mass()
    molecule1.rotate_on_principal_axes()
    molecule2.rotate_on_principal_axes()

    # move the 
    molecule2.translate(distance, 0, 0)

    for attempt in range(attempts):

        molecule1.rotate(*get_random_angle(random))
        molecule2.rotate(*get_random_angle(random), distance, 0, 0)
        # move molecule2 to the minimum distance away from molecule1

        closest_distance = min_inter_distance
        for atom1 in molecule1.get_atoms():
            for atom2 in molecule2.get_atoms():
                if atom1.distance(atom2) < closest_distance:
                    closest_distance = atom1.distance(atom2)

        if not closest_distance < min_inter_distance:
            return

    raise RanOutOfAttemptsException

class RanOutOfAttemptsException(Exception):
    pass


def get_random_angle(random):
    x_angle = math.pi - random.random() * math.pi*2
    y_angle = math.asin(random.random()) * -1 if random.random() < 0.5 else 1
    z_angle = 0
    return x_angle, y_angle, z_angle

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python configuration_generator_2b.py <monomer1xyz> <monomer2xyz> <num configs> <config file>")
        exit(1)

    generate_configurations(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
