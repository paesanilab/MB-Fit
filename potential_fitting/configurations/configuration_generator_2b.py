# External package imports
import math
from random import Random, randint

# Absolute module impots
from potential_fitting.utils import Quaternion
from potential_fitting.molecule import Molecule, xyz_to_molecules
from potential_fitting.utils import files, SettingsReader

def generate_2b_configurations(settings_path, geo1_path, geo2_path, number_of_configs, config_path, min_distance = 1, 
        max_distance = 5, min_inter_distance = 0.8, progression = False, use_grid = False, step_size = 0.5,
        num_attempts = 100, logarithmic = False, seed = None):
        
    """
    Generates a set of 2 body configurations of the two optimized geometries and outputs them to an xyz file.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        geo1_path           - Local path to the first optimized (or series of unoptimized) geometry.
        geo2_path           - Local path to the second optimized (or series of unoptimized) geometry.
        number_of_configs   - Number of configurations to generate.
        config_path         - Local path to the file to write the configurations.
        min_distance        - The minimum distance between the centers of mass of the two monomers.
        max_distance        - The maximum distance between the centers of mass of the two monomers.
        min_inter_distance  - Minimum intermolecular distance is this times the sum of the van der walls radii of two
                atoms.
        progression         - If True, use a linear progression for distance between the centers of masses of the
                two molecules, otherwise use a random progression.
        use_grid            - If True, then distance between the center of mass of the monomers will be on a grid. 
                Otherwise, it will be smooth.
        step_size           - Only used if use_grid is True, this is the step size of the grid in angstroms.
        num_attempts        - The number of attempts to generate a configuration at any given distance before giving
                up and moving to the next distance.
        logarithmic         - If True, then a logarithmic progression is used to generate the configurations.
                This means more configs are generated at lower distances.
        seed                - Seed to use, the same seed will give the same configurations.

    Returns:
        None
    """    
    
    if progression == False:
        generate_2b_configurations_random(settings_path, geo1_path, geo2_path, number_of_configs, config_path, min_distance,
                max_distance, min_inter_distance, num_attempts = num_attempts, logarithmic = logarithmic, seed = seed)
    else:
        generate_2b_configurations_smooth(settings_path, geo1_path, geo2_path, number_of_configs, config_path, min_distance,
                max_distance, min_inter_distance, use_grid = use_grid, step_size = step_size, num_attempts = num_attempts, logarithmic = logarithmic, seed = seed)

def generate_2b_configurations_random(settings_path, geo1_path, geo2_path, number_of_configs, config_path, min_distance = 1, 
        max_distance = 5, min_inter_distance = 0.8, num_attempts = 100, logarithmic = False, seed = None):

    """
    Helper Function to Generate a set of 2 body configurations of the two optimized geometries at random lengths and
    outputs them to an xyz file.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        geo1_path           - Local path to the first optimized (or series of unoptimized) geometry.
        geo2_path           - Local path to the second optimized (or series of unoptimized) geometry.
        number_of_configs   - Number of configurations to generate.
        config_path         - Local path to the file to write the configurations.
        min_distance        - The minimum distance between the centers of mass of the two monomers.
        max_distance        - The maximum distance between the centers of mass of the two monomers.
        min_inter_distance  - Minimum intermolecular distance is this times the sum of the van der walls radii of two
                atoms.
        num_attempts        - The number of attempts to generate a configuration at any given distance before giving
                up and moving to the next distance.
        logarithmic         - If True, then a logarithmic progression is used to generate the configurations.
                This means more configs are generated at lower distances.
        seed                - Seed to use, the same seed will give the same configurations.

    Returns:
        None
    """

    settings = SettingsReader(settings_path)
    
    if seed is None:
        seed = randint(-100000, 100000)

    #parse the molecules from the input xyz files
    molecules1 = xyz_to_molecules(geo1_path)
    molecules2 = xyz_to_molecules(geo2_path)

    import copy

    molecule_init = copy.deepcopy(molecules1[0])

    molecule_init.move_to_center_of_mass()
    molecule_init.rotate_on_principal_axes()

    molecule_rotated = copy.deepcopy(molecules1[0])
    molecule_rotated.move_to_center_of_mass()
    molecule_rotated.rotate_on_principal_axes()
    
    #print("BASE_RMSD:", molecule_init.rmsd2(molecule_rotated))

    # construct a psuedo-random number generator
    random = Random(seed)
    
    #setting the total number of configs
    total_configs = number_of_configs
   
    # open the config file to write to
    with open(files.init_file(config_path, files.OverwriteMethod.get_from_settings(settings)), "w") as config_file:

     
        while total_configs > 0:

            #random distance between min_distance and max_distance

            if logarithmic:
                random_distance = logarithmic_progression(min_distance, max_distance, number_of_configs, random.uniform(0, number_of_configs - 1))
            else:
                random_distance = random.uniform(min_distance, max_distance)
            
            #getting the molecules
            molecule1 = random.choice(molecules1)
            molecule2 = random.choice(molecules2)

            #generating one confiugration at that random distance
        
            try:
                move_to_config(random, molecule1, molecule2, random_distance, min_inter_distance, num_attempts)
            except RanOutOfAttemptsException:
                # if we didn't find a valid configuration, skip this config
                continue

            molecule_rotated = copy.deepcopy(molecule1)
            molecule_rotated.move_to_center_of_mass()
            molecule_rotated.rotate_on_principal_axes()
            
            #print("AFTER RMSD:", molecule_init.rmsd2(molecule_rotated))
            # write total number of atoms to config file
            config_file.write("{}\n".format(molecule1.get_num_atoms() + molecule2.get_num_atoms()))

            # in the comment line, write how many configs have been generated before this one
            config_file.write("{}\n".format(number_of_configs - total_configs))

            # write the xyz of each monomer to the config file
            config_file.write("{}\n{}\n".format(molecule1.to_xyz(), molecule2.to_xyz()))
            
            #decrementing required number of configs
            total_configs -= 1
            
            
    # if we have hit our target number of configs, return
    print("Generated {} configurations".format(number_of_configs))
    return 

def generate_2b_configurations_smooth(settings_path, geo1_path, geo2_path, number_of_configs, config_path, min_distance = 1, 
        max_distance = 5, min_inter_distance = 0.8, use_grid = False, step_size = 0.5, num_attempts = 100,
        logarithmic = False, seed = None):
    """
    Helper Function to Generate a set of 2 body configurations of the two optimized geometries based on a smooth
    progression and outputs them to an xyz file.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        geo1_path           - Local path to the first optimized (or series of unoptimized) geometry.
        geo2_path           - Local path to the second optimized (or series of unoptimized) geometry.
        number_of_configs   - Number of configurations to generate.
        config_path         - Local path to the file to write the configurations.
        min_distance        - The minimum distance between the centers of mass of the two monomers.
        max_distance        - The maximum distance between the centers of mass of the two monomers.
        min_inter_distance  - Minimum intermolecular distance is this times the sum of the van der walls radii of two
                atoms.
        use_grid            - If True, then distance between the center of mass of the monomers will be on a grid. 
                Otherwise, it will be smooth.
        step_size           - Only used if use_grid is True, this is the step size of the grid in angstroms.
        num_attempts        - The number of attempts to generate a configuration at any given distance before giving
                up and moving to the next distance.
        logarithmic         - If True, then a logarithmic progression is used to generate the configurations.
                This means more configs are generated at lower distances.
        seed                - Seed to use, the same seed will give the same configurations.

    Returns:
        None.
    """

    settings = SettingsReader(settings_path)

    if seed is None:
        seed = randint(-100000, 100000)

    # if use_grid is false, set the step size to even space the configurations
    if not use_grid:
        step_size = (max_distance - min_distance) / number_of_configs

    # how many steps the grid will have
    num_steps = math.floor((max_distance - min_distance) / step_size)

    # how many configurations to generate per step in the grid (will be 1 if use_grid is False)
    configs_per_distance = math.ceil(number_of_configs / num_steps)

    # parse the molecules from the input xyz files
    molecules1 = xyz_to_molecules(geo1_path)
    molecules2 = xyz_to_molecules(geo2_path)

    # keeps track of how many total configurations have been generated
    total_configs = 0

    # construct a psuedo-random number generator
    random = Random(seed)
    
    # open the config file to write to
    with open(files.init_file(config_path, files.OverwriteMethod.get_from_settings(settings)), "w") as config_file:
    
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
                    if logarithmic:
                        move_to_config(random, molecule1, molecule2, logarithmic_progression(min_distance, max_distance, number_of_configs, number_of_configs * step / (num_steps)), min_inter_distance, num_attempts)
                    else:
                        move_to_config(random, molecule1, molecule2, min_distance + step * step_size, min_inter_distance, num_attempts)

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

def move_to_config(random, molecule1, molecule2, distance, min_inter_distance, attempts):
    """
    Moves the molecules to a configuration with the given distance between their centers of
    mass. Raises RanOutOfAttemptsException if a configuration is failed to be found after a certain number of attempts.

    Args:
        random              - The Random object used to generate the configuration.
        molecule1           - Molecule object of the 1st monomer.
        molecule2           - Molecule object of the 2nd monomer.
        distance            - Distance between the centers of mass of the two molecules.
        min_inter_distance  - Minimum intermolecular distance is this times the sum of the atoms van der walls radii.
        attempts            - The number of configurations to try before giving up and raising the Exception.

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

    for attempt in range(attempts):

        # rotate each molecule a random amount
        molecule1.rotate(Quaternion.get_random_rotation_quaternion(random), 0, 0, 0)
        molecule2.rotate(Quaternion.get_random_rotation_quaternion(random), distance, 0, 0)

        # setting a flag variable to keep track of a valid configuration
        flag = True


        #loop through atoms in molecule 1
        for atom1 in molecule1.get_atoms():
            #loop through atoms in molecule 2
            for atom2 in molecule2.get_atoms():
                #get the sum of the venderwall radii of the pair
                sum_vdw_distance = atom1.get_vdw_radius() + atom2.get_vdw_radius()
                #comparing with our factor (default = 0.8)
                if atom1.distance(atom2) < min_inter_distance * sum_vdw_distance:
                    flag = False

        #Returning only valid configurations.             
        if flag:
            return

    # if we run out of attempts without generating a valid configuration, raise an exception
    raise RanOutOfAttemptsException
    
class RanOutOfAttemptsException(Exception):
    """
    Used to check if move_to_config runs out of attempts
    """
    pass
