import sys, os
import math
from random import Random

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../qm_mb_energy_calculator/src/")

from molecule import Molecule

def generate_rigid_configurations(geo1, geo2, number_of_configs, config_path, min_distance = 1, max_distance = 10, min_inter_distance = 1.6, use_grid = False, step_size = 0.5):
    """
    Generates a set of 2 body configurations of the two optimized geometries and outputs them to an xyz file

    Args:
        geo1        - the first optimized geometry
        geo2        - the second optimized geometry
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

    num_steps = math.ceil((max_distance - min_distance) / step_size)

    configs_per_distance = math.ceil(number_of_configs / num_steps)

    with open(geo1, "r") as geometry1_file:
        molecule1 = Molecule().read_xyz(geometry1_file, ["CO2"], [3], [0], [1])

    with open(geo2, "r") as geometry2_file:
        molecule2 = Molecule().read_xyz(geometry2_file, ["O3"], [3], [0], [1])

    # put both geoemtries in standard orientation if they are not already
    molecule1.move_to_center_of_mass()
    molecule2.move_to_center_of_mass()
    molecule1.rotate_on_principal_axes()
    molecule2.rotate_on_principal_axes()

    # move molecule2 to the minimum distance away from molecule1
    molecule2.translate(min_distance, 0, 0)

    total_configs = 0

    # construct a psuedo-random number generator
    random = Random()

    print(step_size)

    with open(config_path, "w") as config_file:

        for step in range(num_steps):

            for config in range(configs_per_distance): 

                molecule1.rotate(*get_random_angle(random), 0, 0, 0)
                molecule2.rotate(*get_random_angle(random), min_distance + step * step_size, 0, 0)

                closest_distance = min_inter_distance
                for atom1 in molecule1.get_atoms():
                    for atom2 in molecule2.get_atoms():
                        if atom1.distance(atom2) < closest_distance:
                            closest_distance = atom1.distance(atom2)

                if closest_distance < min_inter_distance:
                    continue

                config_file.write("{}\n".format(molecule1.get_num_atoms() + molecule2.get_num_atoms()))

                config_file.write("{}\n".format(total_configs))

                config_file.write("{}\n{}\n".format(molecule1.to_xyz(), molecule2.to_xyz()))
            
                total_configs += 1

                if total_configs == number_of_configs:
                    return
            # tranlate molecule2 one step further away from molecule1
            molecule2.translate(step_size, 0, 0)

def get_random_angle(random):
    x_angle = math.pi - random.random() * math.pi*2
    y_angle = math.asin(random.random()) * -1 if random.random() < 0.5 else 1
    z_angle = 0
    return x_angle, y_angle, z_angle

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python 2b_configuration_generator.py <monomer1xyz> <monomer2xyz> <num configs> <config file>")
        exit(1)

    generate_rigid_configurations(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
