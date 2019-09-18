from random import Random

from potential_fitting.utils import Quaternion, system
from potential_fitting.molecule import Molecule

from .configuration_generator import ConfigurationGenerator

class AtomDistanceConfigurationGenerator(ConfigurationGenerator):
    """
    Implementation of ConfigurationGenerator that generates configurations using
    a random sampling of rotations while holding two atoms at a certain distance.
    """

    def __init__(self, settings_path, mol1_atom_index, mol2_atom_index, distance=2, min_inter_distance=0.8, num_attempts=100):
        """
        Constructs a new AtomDistanceConfigurationGenerator.

        Args:
            settings_path       - Local path to '.ini' settings file with all relevant settings.
            mol1_atom_index     - The index of the atom in the first monomer to rotate around.
            mol2_atom_index     - The index of the atom in the second monomer to rotate around.
            distance            - The distance between the two atoms.
            min_inter_distance  - Minimum intermolecular distance is this times the sum of the van der walls radii of two
                    atoms.
            num_attempts        - The number of attempts to generate a configuration at any given distance before giving
                    up and moving to the next distance.

        Returns:
            A new AtomDistanceConfigurationGenerator.
        """

        super(AtomDistanceConfigurationGenerator, self).__init__(settings_path)

        self.mol1_atom_index = mol1_atom_index
        self.mol2_atom_index = mol2_atom_index
        self.distance = distance
        self.min_inter_distance = min_inter_distance
        self.num_attempts = num_attempts

    def move_to_config(self, random, molecule1, molecule2, distance):
        """
        Moves the molecules to a configuration with the given distance between their atoms of rotation.

        Raises RanOutOfAttemptsException if a configuration is failed to be found after a certain number of attempts.

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

        mol1_atom = molecule1.get_atoms()[self.mol1_atom_index]
        mol2_atom = molecule2.get_atoms()[self.mol2_atom_index]

        molecule1.translate(-mol1_atom.get_x(), -mol1_atom.get_y(), -mol1_atom.get_z())
        molecule2.translate(-mol2_atom.get_x(), -mol2_atom.get_y(), -mol2_atom.get_z())

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
        Generates Configurations by rotating each molecule around the atom specified in the constructor
        while holding those atoms at a specific distance.

        Might generate less than num_configs configs if min_inter_distance is set to large.

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

        molecules1 = molecule_lists[0]
        molecules2 = molecule_lists[1]

        system.format_print(
            "Beginning 2B configurations generation with constant distance between atoms at indices {} and {}.".format(
                self.mol1_atom_index, self.mol2_atom_index),
            bold=True, color=system.Color.YELLOW)

        total_configs = 0

        for i in range(num_configs):

            molecule1 = random.choice(molecules1)
            molecule2 = random.choice(molecules2)

            try:
                self.move_to_config(random, molecule1, molecule2, self.distance)
            except RanOutOfAttemptsException:
                continue

            mol = Molecule.read_xyz_direct(str(molecule1.get_num_atoms() + molecule2.get_num_atoms()) + "\n\n" + molecule1.to_xyz() + "\n" + molecule2.to_xyz())

            yield mol

            total_configs += 1

            if total_configs % 100 == 0:
                system.format_print("{} configs done...".format(total_configs),
                                    italics=True)

        system.format_print("Done! Generated {} configurations.".format(total_configs), bold=True,
                                color=system.Color.GREEN)

        if total_configs < num_configs:
            system.format_print("Generated fewer than {} configs because it was too hard to generate configurations at the given distance without placing atoms too close together.".format(num_configs), bold=True,
                                color=system.Color.GREEN)

class RanOutOfAttemptsException(Exception):
    """
    Used to check if move_to_config runs out of attempts
    """
    pass