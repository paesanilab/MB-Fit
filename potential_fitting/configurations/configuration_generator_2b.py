import math
from random import Random, randint

from potential_fitting.utils import Quaternion

from potential_fitting.molecule import Molecule, xyz_to_molecules

def generate_2b_configurations_random(geo1, geo2, number_of_configs, config_path, min_distance = 1, 
        max_distance = 5, min_inter_distance = 0.8, seed = randint(-100000, 1000000)):

    #parse the molecules from the input xyz files
    molecules1 = xyz_to_molecules(geo1)
    molecules2 = xyz_to_molecules(geo2)

    # construct a psuedo-random number generator
    random = Random(seed)
    
    #num_attempts is set to five (default)
    num_attempts = 5
    
    #setting the total number of configs
    total_configs = number_of_configs

   
    # open the config file to write to
    with open(config_path, "w") as config_file:

     
        while total_configs > 0:

            #random distance between min_distance and max_distance
            
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



def generate_2b_configurations_smooth(geo1, geo2, number_of_configs, config_path, min_distance = 1, 
        max_distance = 5, min_inter_distance = 0.8, use_grid = False, step_size = 0.5, seed = randint(-100000, 1000000)):
    """
    Generates a set of 2 body configurations of the two optimized geometries and outputs them to an xyz file

    Args:
        geo1        - the first optimized (or series of unoptimized) geometry
        geo2        - the second optimized (or series of unoptimized) geometry
        number_of_configs - generate this many configurations
        config_path - the path to the file to write the configurations
        min_distance - the minimum distance between the centers of mass of the two monomers
        max_distance - the maximum distance between the centers of mass of the two monomers
        min_inter_distance - the factor to multiply the sum of the vanderwall_radius 
                             of the 2 body configuration.
        use_grid    - if True, then distance between the center of mass of the monomers will be on a grid, otherwise, it will be smooth
        step_size - only used if use_grid is True, this is the step size of the grid in angstroms
        num_attempts - the number of attempts before giving up.

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
    
    #num_attempts is set to five (default)
    num_attempts = 5

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




def move_to_config(random, molecule1, molecule2, distance, min_inter_distance, attempts):
    """
    Moves the given molecules to a configuration with the given distance between their centers of
    mass. Raises RanOutOfAttemptsException if a configuration is failed to be found.

    Args:
        random      - the Random object used to generate the configuration
        molecule1   - the Molecule object of the 1st monomer
        molecule2   - the Molecule object of the 2nd monomer
        distance    - distance between the centers of mass of the two molecules
        min_inter_distance - the factor to multiply the sum of the vanderwall_radius 
                             of the 2 body configuration.
        attempts    - the number of configurations to try before giving up

    Returns:
        None

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

        # calculate the minimum distance of any intermolecular interaction
        flag = True



        for atom1 in molecule1.get_atoms():
            for atom2 in molecule2.get_atoms():
                sum_vdw_distance = atom1.get_vdw_radius() + atom2.get_vdw_radius()
                if atom1.distance(atom2) < min_inter_distance * sum_vdw_distance:
                    flag = False

        if flag:
            return 


        
        '''Ethan's Old Code, will change once final changes are done.

        #closest_distance = min_inter_distance
        #for atom1 in molecule1.get_atoms():
            #for atom2 in molecule2.get_atoms():
                #if atom1.distance(atom2) < closest_distance:
                    #closest_distance = atom1.distance(atom2)

        # if the minimum intermolecular distance in this configuration is ngreater than or equal to min_inter_distance, then this is a valid configuration

        #if closest_distance >= min_inter_distance:
            #return

        # otherwise, we repeat the loop and make another attempt
        '''

    # if we run out of attempts without generating a valid configuration, raise an exception
    raise RanOutOfAttemptsException
    



def generate_2b_configurations(geo1, geo2, number_of_configs, config_path, min_distance = 1, 
        max_distance = 5, min_inter_distance = 0.8, progression = False, use_grid = False, step_size = 0.5, seed = randint(-100000, 1000000)):

    if progression == False:
        generate_2b_configurations_random(geo1, geo2, number_of_configs, config_path, min_distance, max_distance, min_inter_distance, seed)
    else:
        generate_2b_configurations_smooth(geo1, geo2, number_of_configs, config_path, min_distance, max_distance, min_inter_distance, use_grid = False, step_size = 0.5, seed = randint(-100000, 1000000))


     



class RanOutOfAttemptsException(Exception):
    """
    Used to check if move_to_config runs out of attempts
    """
    pass
