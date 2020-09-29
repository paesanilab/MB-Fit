import math, numpy
from random import Random

from mbfit.utils import system, constants
from mbfit.exceptions import LineFormatError, ParsingError, InvalidValueError, InconsistentValueError
from mbfit.utils.distribution_function import ConstantDistributionFunction, PiecewiseDistributionFunction,\
                                                           LinearDistributionFunction, GeometricDistributionFunction

from .configuration_generator import ConfigurationGenerator



class NormalModesConfigurationGenerator(ConfigurationGenerator):
    """
    Implementation of ConfigurationGenerator that generates configurations from normal
    mode data.
    """

    def __init__(self, settings_path, normal_modes_path, classical=True, distribution = 'piecewise',
        temperature=None, distribution_function=None):
        """
        Constructs a new NormalModesConfigurationGenerator.

        Args:
            settings_path       - Local path to '.ini' settings file with all relevant settings.
            normal_modes_path   - Local path to the '.dat' file containing normal modes information.
            classical           - If True, use a classical distribution to generate the configurations, otherwise quantum.
            distribution        - One of the following choices: 'piecewise', 'constant', 'linear', 'geometric', 'custom'
                    'piecewise' uses a piecewise distribution in the following style:
                        5% at highest frequency / 100
                        40% at highest frequency / 20
                        30% at highest frequency / 10
                        20% at highest frequency / 5
                        5% at highest frequency / 2
                    'constant' uses a set temperature for all configurations.
                        Specify the temperature by setting the temperature argument to a single value.
                    'linear' uses a linear distribution from a minimum to a maximum temperature.
                        Specify the min and max temperature by setting the temperature argument to a 2-tuple: (min, max).
                        If temperature is unspecified, then the minimum is 0, and the maximum is the highest normal mode frequency.
                    'geometric' uses a geometric distribution from a minimum to a maximum temperature.
                        Specify the min and max temperature by setting the temperature argument to a 2-tuple: (min, max).
                        If temperature is unspecified, then the minimum is 0, and the maximum is the highest normal mode frequency.
                    'custom' uses a user-specified DistributionFunction to generate the temperatures used during configuration generation.
            temperature         - Should be set to different values based on what distribution is being used.
                    If 'piecewise' or 'custom' distribution, then temperature is ignored.
                    If 'constant' distribution, then temperature should be a single value.
                    If 'linear' or 'geometric' distribution, then temperature should be a 2-tuple: (min, max)
                    All temperatures should be specified in KELVIN.
            distribution_function - Implementation of DistributionFunction. Only used if distribution='custom'.
                    distribution_function.get_vale(x) should be implemented over the domain [0,1]. So the first config
                    will have temperature = distribution_function.get_value(0) and the last config will have temperature =
                    distribution_function.get_value(1), with configurations in between passing linearly increasing values to 
                    distribution_function.get_value(x).
                    The distribution_function should return temperatures in atomic units (NOT KELVIN).
                    See package utils.distribution_function for abstract DistributionFunction class and example implementaitons.

        Returns:
            A new NormalModesConfigurationGenerator.
        """

        super(NormalModesConfigurationGenerator, self).__init__(settings_path)

        # parse frequencies, reduced masses, and normal modes from the input file.
        self.frequencies, self.reduced_masses, self.normal_modes = self.parse_normal_modes_file(normal_modes_path)

        # sort frequencies, reduced masses, and normal modes from smallest frequency to largest.
        self.frequencies, self.reduced_masses, self.normal_modes = self.sort_by_frequency(self.frequencies, self.reduced_masses, self.normal_modes)

        # Check for negative frequencies and print warnings.
        num_neg_freqs = 0
        for frequency in self.frequencies:
            if frequency < 0:
                num_neg_freqs += 1

        if num_neg_freqs == 1:
            system.format_print(
                "Single negative frequency detected in input. This most likely means the given geometry is a transition state.",
                italics=True)

        elif num_neg_freqs > 1:
            system.format_print(
                "Multiple ({}) negative frequencies detected in input. Proceed with caution.".format(num_neg_freqs),
                italics=True)

        # convert the frequencies to atomic units from cm-1
        # absolute value allows supporting of imaginary frequencies (i think)
        self.frequencies = [abs(frequency) / constants.autocm for frequency in self.frequencies]

        if distribution == 'piecewise':

            system.format_print("Will generate configurations over a piecewise temperature distribution.",
                                italics=True)
            self.temp_distribution = PiecewiseDistributionFunction(
                [
                ConstantDistributionFunction(self.frequencies[-1] / 100),
                ConstantDistributionFunction(self.frequencies[-1] / 20),
                ConstantDistributionFunction(self.frequencies[-1] / 10),
                ConstantDistributionFunction(self.frequencies[-1] / 5),
                ConstantDistributionFunction(self.frequencies[-1] / 2),
                ],
                [
                0.05, 0.45, 0.75, 0.95
                ]
            )
        elif distribution == 'constant':

            if temperature == None:
                raise InvalidValueError('temperature', temperature, "please specify a temperature for the constant distribution")

            system.format_print("Will generate configurations at temperature {} K.".format(temperature),
                                italics=True)

            temperature *= constants.kelvin_to_au
            self.temp_distribution = ConstantDistributionFunction(temperature)

        elif distribution == 'linear':


            if temperature == None:
                temperature = [0, self.frequencies[-1] / constants.kelvin_to_au]

            system.format_print("Will generate configurations over a linear temperature distribution from {} K to {} K.".format(temperature[0], temperature[1]),
                                italics=True)
            self.temp_distribution = LinearDistributionFunction.get_function_from_2_points(0, temperature[0] * constants.kelvin_to_au,
                                                                                           1, temperature[1] * constants.kelvin_to_au)

        elif distribution == 'geometric':

            if temperature == None:
                temperature = [self.frequencies[0] / constants.kelvin_to_au, self.frequencies[0] * (self.frequencies[0] / (2*self.frequencies[-1])) ** -1 / constants.kelvin_to_au]

            system.format_print("Will generate configurations over a geometric temp distribution from {} K to {} K.".format(temperature[0], temperature[1]),
                                italics=True)
            self.temp_distribution = GeometricDistributionFunction(temperature[0]*constants.kelvin_to_au,
                                                                  temperature[1] / temperature[0])

        elif distribution == 'custom':
            system.format_print("Will generate configurations over a user-specified distribution.",
                                italics=True)

            if distribution_function is None:
                raise InvalidValueError("distribution_function", distribution_function, "because distribution='custom', you must specify a distribution_funciton")

            self.temp_distribution = distribution_function

        else:
            raise InvalidValueError("distribution", distribution, "use one of 'piecewise', 'constant', 'linear', 'geometric', or 'custom'")

        self.classical = classical

        system.format_print("Temp Distribution: {} for x in range [0,1].".format(self.temp_distribution.to_string(dep_name="temp (au)")),
                        italics=True)

    def parse_normal_modes_file(self, normal_modes_path):
        """
        Reads a normal modes file and parses the frequencies, reduced masses, and normal modes from it.

        Args:
            normal_modes_path   - Local path to the '.dat' file containing normal modes information.

        Returns:
            (frequencies, reduced_masses, normal_modes) as parsed from input file.
        """

        system.format_print("Parsing normal modes input file {}".format(normal_modes_path), italics=True)

        frequencies = []
        reduced_masses = []
        normal_modes = []

        num_atoms = sum([int(atom_num) for atom_num in self.settings.get("molecule", "fragments").split(",")])

        # read frequencies, reduced masses, and normal modes from the input file.
        with open(normal_modes_path, "r") as normal_modes_file:

            # we loop until we run out of normal modes to parse
            while (True):
                first_line = normal_modes_file.readline()

                # if first line is empty string, we have reached EOF
                if first_line == "":
                    break

                # if first line is not of valid format, raise error
                if not first_line.startswith("normal mode:") or not len(first_line.split()) == 3:
                    raise LineFormatError(normal_modes_path, first_line, "EOF or normal mode: x")

                frequency_line = normal_modes_file.readline()

                # if the frequency line is not of valid format, raise error
                if frequency_line == "":
                    raise ParsingError(normal_modes_path, "Unexpected EOF, expected line of format 'frequency = x'")

                if not frequency_line.startswith("frequency = ") or not len(frequency_line.split()) == 3:
                    raise LineFormatError(normal_modes_path, frequency_line, "frequency = x")

                # parse the frequency from the frequency line
                try:
                    frequency = float(frequency_line.split()[2])

                except ValueError:
                    raise ParsingError(normal_modes_path,
                                       "cannot parse {} into a frequency float".format(
                                           frequency_line.split()[2])) from None
                frequencies.append(frequency)

                reduced_mass_line = normal_modes_file.readline()

                # if the reduced mass line is not of valid format, raise error
                if reduced_mass_line == "":
                    raise ParsingError(normal_modes_path, "Unexpected EOF, expected line of format 'recued mass = x'")

                if not reduced_mass_line.startswith("reduced mass = ") or not len(reduced_mass_line.split()) == 4:
                    raise LineFormatError(normal_modes_path, reduced_mass_line, "reduced mass = x")

                # parse the reduced mass from the frequency line
                try:
                    reduced_mass = float(reduced_mass_line.split()[3])

                except ValueError:
                    raise ParsingError(normal_modes_path,
                                       "cannot parse {} into a reduced mass float".format(
                                           reduced_mass_line.split()[3])) from None

                reduced_masses.append(reduced_mass)

                # parse the normal mode
                normal_mode = [[None, None, None] for i in range(num_atoms)]

                for atom_index in range(num_atoms):
                    normal_mode_line = normal_modes_file.readline()

                    if normal_mode_line == "":
                        raise ParsingError(normal_modes_path,
                                           "Unexpected EOF, expected line of format 'x y z'")

                    if len(normal_mode_line.split()) != 3:
                        raise LineFormatError(normal_modes_path, normal_mode_line,
                                              "x y z")

                    for ordinate_index, token in enumerate(normal_mode_line.split()):
                        try:
                            offset = float(token)

                        except ValueError:
                            raise ParsingError(normal_modes_path,
                                               "cannot parse {} into a offset float".format(token)) from None

                        normal_mode[atom_index][ordinate_index] = offset

                normal_modes.append(normal_mode)

                # skip the blank line
                blank_line = normal_modes_file.readline()
                if blank_line != "\n":
                    raise ParsingError(normal_modes_path, "expected blank line")

        system.format_print("Completed parsing normal modes input file.", italics=True)

        return frequencies, reduced_masses, normal_modes

    def sort_by_frequency(self, frequencies, reduced_masses, normal_modes):
        """
        Sorts the given lists of frequencies, reduced masses, and normal modes by their frequencies from least to greatest.

        The input lists must be ordered such that normal mode x has frequency frequencies[x], reduced mass
        reduced_masses[x] and offsets normal_modes[x].

        Args:
            frequencies         - Frequencies of the normal modes.
            reduced_masses      - Reduced masses of the normal modes.
            normal_modes        - Offset of each atom for each normal mode. Each element of this list contains one sub-list
                    for each atom in the molecule. Each of these sublists is [x offset, y offset, z offset] of the normal
                    mode.

        Returns:
            A 3-tuple (frequencies, reduced_masses, normal_modes) sorted by frequencies.
        """

        reduced_masses = [mass for freq, mass in sorted(zip(frequencies, reduced_masses))]
        normal_modes = [normal_mode for freq, normal_mode in sorted(zip(frequencies, normal_modes))]
        frequencies = sorted(frequencies)

        return frequencies, reduced_masses, normal_modes

    def make_config(self, molecule, G, random):
        """
        Gets a single configuration based on the input molecule, and G.

        Args:
            molecule            - The molecule to generate a configuration of.
            G                   - The sqrt of the mass-scaled covariance matrix.
            random              - The random object to use to generate the configuration.

        Returns:
            A configurations
        """

        molecule = molecule.get_copy()

        # calculate the dimension of this molecule
        dim = 3 * molecule.get_num_atoms()

        # initialize the displacement lists to all zero, they track this configurations displacement from the
        # optimized geometry
        displacement = [[0, 0, 0] for i in range(molecule.get_num_atoms())]

        # generate a list of random numbers in a normal distribution, with mean 0 and standard deviation 1
        norm_dist_list = [random.normalvariate(0, 1) for i in range(dim)]

        # loop over each atom's displacement
        for atom_index, atom, atom_displacement in zip(range(molecule.get_num_atoms()), molecule.get_atoms(),
                                                       displacement):

            # loop over x, y, and z in the current atom's displacement
            for coordinate_index in range(3):
                # set the displacement equal to the inner product of the corresponding column vector from G and the
                # random number list
                atom_displacement[coordinate_index] = numpy.dot([g[atom_index * 3 + coordinate_index] for g in G],
                                                                norm_dist_list)

                # unscale the atom displacement ordinate by the molecules mass relative to that of an electron
                atom_displacement[coordinate_index] /= math.sqrt(
                    atom.get_mass() * constants.mass_electron_per_mass_proton)

        # scale the bohr constants from meters to angstroms
        bohr = constants.bohr * 1e10

        # loop over each atom in the molecule
        for atom_index, atom in enumerate(molecule.get_atoms()):
            # scale each atom's coordinates to atomic units and add the displacement
            x = atom.get_x() / bohr + displacement[atom_index][0]
            y = atom.get_y() / bohr + displacement[atom_index][1]
            z = atom.get_z() / bohr + displacement[atom_index][2]

            # convert back to angstroms
            x *= bohr
            y *= bohr
            z *= bohr

            atom.set_xyz(x, y, z)

        return molecule

    def generate_configurations(self, molecule_lists, num_configs, seed=None):
        """
        Generates Normal modes configurations of the given molecule.

        Args:
            molecule_lists  - A List of lists containing only a single element such that molecule_lists[0][0] is the
                    optimized geometry for the configuration generation.
            num_configs     - The number of configurations to generate.
            seed            - Seed for the random number generator. The same seed will yield the same configurations
                    when all else is held equal.

        Yields:
            Molecule objects containing the new configurations.
        """

        system.format_print("Beginning normal modes configuration generation.",
                            bold=True, color=system.Color.YELLOW)

        if seed is None:
            seed = self.get_rand_seed()

        # parse the molecule from the input ".xyz" into a Molecule object
        molecule = molecule_lists[0][0]

        # create a new random object from the seed.
        random = Random(seed)

        # calculate the dimension of this molecule
        dim = 3 * molecule.get_num_atoms()
        # calculate the dimension of the null space of this molecule
        dim_null = dim - len(self.normal_modes)

        # check to make sure the number of normal modes is right.
        if dim_null < 5:
            raise InvalidValueError("number of normal modes", dim - dim_null,
                                    "larger than 3 * (number of atoms in the molecule) - 5 ({})".format(dim))

        system.format_print("Will generate {} configs over the temperature distribution.".format(num_configs),
                            italics=True)

        # mass-scale and normalize the normal modes
        for normal_mode in self.normal_modes:

            # normalization scale is used to keep track of the length of the normal mode vector to normalize it.
            normalization_scale = 0

            for coordinates, atom in zip(normal_mode, molecule.get_atoms()):

                # scale each element of the normal modes relative to its atom's mass relative to the mass of an electron
                # (for some reason)
                sqrt_mass = math.sqrt(atom.get_mass() * constants.mass_electron_per_mass_proton)
                for i in range(3):
                    coordinates[i] = coordinates[i] * sqrt_mass

                # add each ordinate squared to the normalization scale
                normalization_scale += coordinates[0] ** 2
                normalization_scale += coordinates[1] ** 2
                normalization_scale += coordinates[2] ** 2

            normalization_scale = math.sqrt(normalization_scale)

            # normalize the normal mode by dividing by the normalization scale
            for coordinates in normal_mode:
                for i in range(3):
                    coordinates[i] = coordinates[i] / normalization_scale

        freq_cutoff = 10 * constants.cmtoau

        # Generate the temp configs.
        if num_configs > 0:
            system.format_print("Generating Temperature Distribution Configs...", italics=True)

            # loop over each temp distribution config to generate
            for config_index in range(num_configs):

                temp = self.temp_distribution.get_value(config_index / (num_configs - 1))

                # fill G with all 0s
                G = [[0 for i in range(dim)] for k in range(dim)]  # sqrt of the mass-scaled covariance matrix

                # for each normal mode, frequency pair, update d and G.
                for normal_mode_index, frequency, reduced_mass, normal_mode in zip(range(len(self.frequencies)),
                                                                                   self.frequencies,
                                                                                   self.reduced_masses,
                                                                                   self.normal_modes):

                    # check if frequency is high enough to have an effect
                    if frequency >= freq_cutoff:

                        if self.classical:
                            d = temp / (frequency ** 2)

                        elif temp > 1.0e-8:
                            # if temp is significantly larger than 0, then set this normal mode's d by the formula
                            d = 0.5 / (numpy.tanh(frequency / (2 * temp)) * frequency)

                        # if temp is not significantly larger than 0 (so it is close to 0), then we must use a different
                        # formula to avoid divide-by-zero error.
                        else:
                            d = 0.5 / frequency

                        # G = ( d * U * U^T )^(1/2), where U are normal modes
                        for i in range(dim):
                            for j in range(dim):
                                G[i][j] += math.sqrt(d) * normal_mode[i // 3][i % 3] * normal_mode[j // 3][j % 3]

                yield self.make_config(molecule, G, random)

            system.format_print("... Successfully generated temperature distribution configs!", italics=True)

#        # Generate the A configs.
#        if num_A_configs > 0:
#
#            system.format_print("Generating A Distribution Configs...", italics=True)
#
#            # loop over each A distribution config to generate
#            for config_index in range(num_A_configs):
#
#               A = self.A_distribution.get_value(config_index / (num_temp_configs - 1))
#
#                # fill G with all 0s
#                G = [[0 for i in range(dim)] for k in range(dim)]  # sqrt of the mass scaled covariance matrix
#
#                # for each normal mode, frequency pair, update d and G.
#                for normal_mode_index, frequency, reduced_mass, normal_mode in zip(range(len(self.frequencies)),
#                                                                                   self.frequencies,
#                                                                                   self.reduced_masses,
#                                                                                   self.normal_modes):
#
#                    # check if frequency is high enough to have an effect
#                    if frequency >= freq_cutoff:
#
#                        if self.classical:
#                            d = temp / (frequency ** 2)
#
#                        # if A is significantly larger than 0, then set this normal mode's d by the formula
#                        if A > 1.0e-8:
#                            d = 0.5 / (numpy.tanh(0.5 / A) * frequency)
#
#                        # if A is not significantly larger than 0 (so it is close to 0), then we must use a different
#                        # formula to avoid divide-by-zero error.
#                        else:
#                            d = 0.5 / frequency
#
#                        # G = ( d * U * U^T )^(1/2), where U are normal modes
#                        for i in range(dim):
#                            for j in range(dim):
#                                G[i][j] += math.sqrt(d) * normal_mode[i // 3][i % 3] * normal_mode[j // 3][
#                                    j % 3]
#
#                yield self.make_config(molecule, G, random)
#
#            system.format_print("... Successfully generated A distribution configs!", italics=True)

        system.format_print("Normal Distribution Configuration generation complete! Generated {} configs.".format(
            num_configs), bold=True, color=system.Color.GREEN)
