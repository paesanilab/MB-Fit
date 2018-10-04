from potential_fitting.molecule import xyz_to_molecules
from potential_fitting.utils import SettingsReader, utils
from potential_fitting.exceptions import ParsingError, LineFormatError, InvalidValueError
from random import randint, Random
from potential_fitting.utils import constants
import math, numpy, copy
def generate_1b_configurations(settings_path, geo_path, normal_modes_path, config_path,
        seed = randint(-1000000, 1000000)):
    """
    Generates a set of 1b configurations for a molecule given its optimized geometry and normal modes.

    Args:
        settings_path - path to the .ini file with relevant settings
        geo_path - path to the optimized geometry .xyz file
        normal_modes_path - path to normal modes file as generated by generate_normal_modes
        config_path - path to the file to write the configurations, existing content will be clobbered.
        seed - use this seed to generate random numbers, the same seed will give the same configurations.
    """

    print("Parsing normal mode input file.")

    settings = SettingsReader(settings_path)

    # read the optimized geometry    
    molecule = xyz_to_molecules(geo_path, settings)[0]

    frequencies = []
    reduced_masses = []
    normal_modes = []

    # read frequencies, reduced masses, and normal modes from the input file.
    with open(normal_modes_path, "r") as normal_modes_file:

        # we loop until we run out of normal modes to parse
        while(True):
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
                        "cannot parse {} into a frequency float".format(frequecy_line.split()[2])) from None
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
                        "cannot parse {} into a reduced mass float".format(reduced_mass_line.split()[3])) from None

            reduced_masses.append(reduced_mass)

            # parse the normal mode
            normal_mode = [[None, None, None] for i in range(molecule.get_num_atoms())]

            for atom_index in range(molecule.get_num_atoms()):
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

                    except ValueErrorL:
                        raise ParsingError(normal_modes_path,
                                "cannot parse {} into a offset float".format(token)) from None

                    normal_mode[atom_index][ordinate_index] = offset
                
            normal_modes.append(normal_mode)

            # skip the blank line
            blank_line = normal_modes_file.readline()
            if blank_line != "\n":
                raise ParsingError(normal_modes_path, "expected blank line")

    print("Completed parsing normal modes input file.")

    generate_1b_normal_mode_configs(settings_path, geo_path, frequencies, reduced_masses, normal_modes, config_path,
            seed = seed)

def generate_1b_normal_mode_configs(settings_path, geo_path, frequencies, reduced_masses, normal_modes, config_path,
        seed = randint(-100000, 100000)):
    """
    Generates 1b configurations based on a geometry (should usually be optimized) and the normal modes corresponding
    to that geometry.

    The first half of the configurations are generated by sampling the harmonic distribution. This distribution
    is over temperature.

    The second half of the configurations are generated by sampling a varation of the harmonic distribution where each
    normal mode is excited based on its frequency. Each mode's temperature is set to A*h_bar*freq/boltzman. This
    distribution is over A.

    If there are an odd number of configs, there will be one more temp config than A config.

    Args:
        settings_path       - Local path to the ".ini" file with all relevant settings.
        geo_path            - Local path to the ".xyz" geoemtry file to read the input geometry from.
        frequencies         - A list of frequencies corresponding to the normal modes, should be ordered from least
                to greatest.
        reduced_masses      - A list of the reduced masses of the normal modes in the same order as the frequencies
                list.
        normal_modes        - A list of the coordiantes of each normal mode in the same order as the frequencies list.
                Each list consists of one length 3 list for each atom in the molecule contianing x, y, and z.
        config_path         - Local path to the ".xyz" geometry file to write the configs to.
        seed                - Seed to use to generate the configurations, the same seed will give the same
                configurations.
    """

    print("Running normal distribution configuration generator...")

    # parse the ".ini" file into a SettingsReader object
    settings = SettingsReader(settings_path)

    # parse the molecule from the input ".xyz" into a Molecule object
    molecule = xyz_to_molecules(geo_path, settings)[0]

    # generate a new random seed
    random = Random(seed)

    # get the number of configurations to generate from settings
    num_configs = settings.getint("config_generator", "num_configs")

    geometric = settings.getboolean("config_generator", "geometric")

    # calculate the dimension of this molecule
    dim = 3 * molecule.get_num_atoms()
    # calculate the dimension of the null space of this molecule
    dim_null = dim - len(normal_modes)

    if dim_null < 0:
        raise InvalidValueError("number of normal modes", dim - dim_null,
                "less than 3 * number of atoms in the molecule ({})".format(dim))

    # number of configs using a distribution over A
    num_A_configs = num_configs // 2
    # number of configs using a distribution over temp
    num_temp_configs = num_configs - num_A_configs

    # deep copy the frequencies, reduced masses, and normal modes input array before we change it, and sort them so
    # that they are all sorted from lowest frequency to highest.
    frequencies, reduced_masses, normal_modes = sort_by_frequency(frequencies, reduced_masses, normal_modes)

    # mass-scale and normalize the normal modes
    for normal_mode in normal_modes:

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
            
    # convert the frequencies to atomic units from cm
    # absolute value allows supporting of imaginary frequencies (i think)
    frequencies = [abs(frequency) / constants.autocm for frequency in frequencies]

    # first we will generate the temp distribution configs

    # if a geometric progression was requested, set min, max, factor, and addend of temp to correspond
    if geometric:

        temp_min = frequencies[0] # Kelvin
        temp_max = 2 * frequencies[-1] # Kelvin
        temp_factor = (temp_min / temp_max) ** (-1 / (num_temp_configs - 1))
        temp_addend = 0

    # if a linear progression was requested, set min, max, factor, and addend of temp to correspond
    else:

        temp_min = 0 # Kelvin
        temp_max = frequencies[-1] # Kelvin
        temp_factor = 1
        temp_addend = (temp_max - temp_min) / (num_temp_configs - 1)

    # initialize temp to the temp minimum, it will be increased each iteration of the loop
    temp = temp_min

    freq_cutoff = 10 * constants.cmtoau

    # open the config file to write configurations to.
    with open(config_path, "w") as config_file:

        # loop over each temp distribution config to generate
        for config_index in range(num_temp_configs):

            # fill G with all 0s
            G = [[0 for i in range(dim)] for k in range(dim)] # sqrt of the mass-scaled covariance matrix

            # for each normal mode, frequency pair, update d and G.
            for normal_mode_index, frequency, reduced_mass, normal_mode in zip(range(len(frequencies)), frequencies,
                    reduced_masses, normal_modes):

                # check if frequency is high enough to have an effect
                if frequency >= freq_cutoff:

                    # if temp is significantly larger than 0, then set this normal mode's d by the formula
                    if temp > 1.0e-8:
                        d = 0.5 / (numpy.tanh(frequency / (2 * temp)) * frequency)

                    # if temp is not significantly larger than 0 (so it is close to 0), then we must use a different
                    # formula to avoid divide-by-zero error.
                    else:
                        d = 0.5 / frequency

                    for i in range(dim):
                        for j in range(dim):
                            G[i][j] += math.sqrt(d) * normal_mode[i // 3][i % 3] * normal_mode[j // 3][j % 3]

            make_config(config_file, config_index + 1, molecule, G, random)

            # increase temp
            temp = temp * temp_factor + temp_addend


    # now we will generate the A distribution configs

    # if a geometric progression was requested, set min, max, factor, and addend of A to correspond
    if geometric:
        A_min = 1
        A_max = 2
        A_factor = (A_min / A_max) ** (-1 / (num_A_configs - 1))
        A_addend = 0

    # if a linear progression was requested, set min, max, factor, and addend of A to correspond
    else:
        A_min = 0
        A_max = 2
        A_factor = 1
        A_addend = (A_max - A_min) / (num_A_configs - 1)

    # initialize A to the A minimum, it will be increased each iteration of the loop
    A = A_min

    # open the config file to write configurations to. Open in append mode so as not to overwrite temp configs.
    with open(config_path, "a") as config_file:

        # loop over each A distribution config to generate
        for config_index in range(num_A_configs):

            # fill G with all 0s
            G = [[0 for i in range(dim)] for k in range(dim)] # ???

            # for each normal mode, frequency pair, update d and G.
            for normal_mode_index, frequency, reduced_mass, normal_mode in zip(range(len(frequencies)), frequencies,
                    reduced_masses, normal_modes):

                # check if frequency is high enough to have an effect
                if frequency >= freq_cutoff:

                    # if A is significantly larger than 0, then set this normal mode's d by the formula
                    if A > 1.0e-8:
                        d = 0.5 / (numpy.tanh(0.5 / A) * frequency)

                    # if A is not significantly larger than 0 (so it is close to 0), then we must use a different
                    # formula to avoid divide-by-zero error.
                    else:
                        d = 0.5 / frequency

                    for i in range(dim):
                        for j in range(dim):
                            G[i][j] += math.sqrt(d) * normal_mode[i // 3][i % 3] * normal_mode[j // 3][j % 3]

            make_config(config_file, config_index + num_temp_configs + 1, molecule, G, random)

            # increase A
            A = A * A_factor + A_addend

    print("Normal Distribution Configuration generation complete.")

def sort_by_frequency(frequencies, reduced_masses, normal_modes):
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

def make_config(config_file, config_index, molecule, G, random):
    """
    Writes a single configuration to the file based on the input molecule, and G.

    Args:
        config_file         - File handle to write the configuration to.
        config_index        - The index of this config in the file. Written in the comment line of the xyz.
        molecule            - The molecule to generate a configuration of.
        G                   - The sqrt of the mass-scaled covariance matrix.
        random              - The random object to use to generate the configuration.

    Returns:
        None.
    """

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

            # de-scale the atom displacement ordinate by the molecules mass relative to that of an electron
            atom_displacement[coordinate_index] /= math.sqrt(atom.get_mass() * constants.mass_electron_per_mass_proton)

    # scale the bohr constants from meters to angstroms
    bohr = constants.bohr*1e10

    # write number of atoms
    config_file.write("{}\n".format(molecule.get_num_atoms()))

    # write index of this config in comment line
    config_file.write("{}\n".format(config_index))

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

        # write this atom's atomic symbol and coordinates
        #config_file.write("{:2} {:22.14e} {:22.14e} {:22.14e}\n".format(atom.get_name(), x, y, z))
        config_file.write("{:2}{:13.8f}{:13.8f}{:13.8f}\n".format(atom.get_name(), x, y, z))
