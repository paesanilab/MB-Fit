# external package imports
import os, sys, contextlib, glob

# local module imports
from .utils import SettingsReader, files, system
from . import configurations, database, polynomials, fitting
from .molecule import xyz_to_molecules
from mbfit.exceptions import InconsistentValueError


def apply_standard_order(settings_path, geo_path):
    """

    Reorders the fragments and atoms in a molecule into 'standard order'.
    This order is arbitrarily defined as:
     * Fragments are ordered alphabetically
     * Atoms within fragments are reordered in the following way:
        * If two atoms have a different atomic number, the atom with the lower number comes first.
        * If two atoms have the same atomic number, the atom with the adjacent atom of the lowest priority comes first.


    WARNING:
    This function is mostly uninteresting and probably not useful unless you know specifically why you are using it.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        geo_path            - Local path to the file containing the non-standard-order xyz data.

    Returns:
        (standard_settings_path, standard_geo_path)
        standard_settings_path  - Local path to new settings file containing information for the molecule in standard order.
            Will be equal to settings_path + "_standard"
        standard_geo_path       - Local path to new xyz file containing data for the molecule in standard order.
            Will be equal to geo_path  + "_standard"
    """

    settings = SettingsReader(settings_path)
    molecule = xyz_to_molecules(geo_path, settings)[0]
    names, fragments, charges, spins, symmetry, SMILES = molecule.get_config_molecule_section()
    settings.set("molecule", "names", names)
    settings.set("molecule", "fragments", fragments)
    settings.set("molecule", "charges", charges)
    settings.set("molecule", "spins", spins)
    settings.set("molecule", "symmetry", symmetry)
    settings.set("molecule", "SMILES", SMILES)

    standard_settings_path = settings_path + "_standard"
    standard_geo_path = geo_path + "_standard"
    
    files.init_file(standard_settings_path, files.OverwriteMethod.BACKUP)
    settings.write(standard_settings_path)

    with open(geo_path, "r") as geo_file:
        geo_file.readline()
        comment_line = geo_file.readline()

    files.init_file(standard_geo_path, files.OverwriteMethod.BACKUP)

    with open(standard_geo_path, "w") as geo_file:
        geo_file.write("{}\n".format(molecule.get_num_atoms()))
        geo_file.write(comment_line)
        geo_file.write("{}\n".format(molecule.to_standard_xyz()))

    return standard_settings_path, standard_geo_path


def optimize_geometry(settings_path, unopt_geo_path, opt_geo_path, method, basis, qm_options={}):
    """
    Optimizes the geometry of the given molecule.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        unopt_geo_path      - Local path to the file to read the unoptimized geoemtry from.
        opt_geo_path        - Local path to the file to write the optimized geometry to.
        method              - The method to use for this geometry optimization.
        basis               - The basis to use for this geometry optimization.
        qm_options           - Dictionary of extra arguments to be passed to the QM code doing the calculation.

    Returns:
        None.
    """

    configurations.optimize_geometry(settings_path, unopt_geo_path, opt_geo_path, method, basis, qm_options=qm_options)

def generate_normal_modes(settings_path, opt_geo_path, normal_modes_path, method, basis, qm_options={}):
    """
    Generates the normal modes for the given molecule.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        opt_geo_path        - Local path to the file to read the optimized geometry from.
        normal_modes_path   - Local path to the file to write the normal modes to.
        method              - The method to use for this normal modes calculation.
        basis               - The basis to use for this normal modes calculation.
        qm_options           - Dictionary of extra arguments to be passed to the QM code doing the calculation.

    Returns:
        Null dimension of normal modes.
    """
    
    dim_null = configurations.generate_normal_modes(settings_path, opt_geo_path, normal_modes_path, method, basis, qm_options=qm_options)

    return dim_null

def generate_normal_mode_configurations(settings_path, opt_geo_path, normal_modes_path, configurations_path,
        number_of_configs=100, seed=None, classical=True, distribution='piecewise',
        temperature=None, distribution_function=None):
    """
    Generates normal mode configurations for a given molecule from a set of normal modes.

    If both linear and geometric are False, will use a piecewise distribution over temperature.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        opt_geo_path        - Local path to the file to read optimized geometry from.
        normal_modes_path   - Local path to the file to read normal modes from.
        dim_null            - The null dimension of this molecule, see generate_normal_modes().
        config_path         - Local path to the file to write configurations to.
        number_of_configs   - Number of configurations to generate
        seed                - The same seed with the same molecule and normal modes will always generate the same
                configurations.
        classical           - If True, use a classical distribution over temp and A, otherwise, use a quantum
                distribution. QM distributions generate a wider distribution over energy.
                Default: True
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
        None.
    """

    config_generator = configurations.NormalModesConfigurationGenerator(settings_path, normal_modes_path,
                                                                        classical=classical,
                                                                        distribution=distribution,
                                                                        temperature=temperature,
                                                                        distribution_function=distribution_function)

    configurations.ConfigurationGenerator.generate_configs_from_file_to_file([opt_geo_path],
                                                                             configurations_path,
                                                                             config_generator,
                                                                             number_of_configs,
                                                                             seed=seed)

def generate_2b_configurations(settings_path, geo1_path, geo2_path, number_of_configs, configurations_path, 
        min_distance=1, max_distance=5, min_inter_distance=0.8, progression=False, use_grid=False,
        step_size=0.5, num_attempts=100, logarithmic=False, distribution=None,
        mol1_atom_index=None, mol2_atom_index=None, seed=None):
    """
    Generates 2b configurations for a given dimer by rotating them randomly over a distribution of
    distances.

    If one of mol1_atom_index is specified, then both must be specified.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        geo1_path           - Local path to read the first optimized geometry from.
        geo2_path           - Local path to read the second optimized geometry from.
        number_of_configs   - The target number of configurations to generate. If max_distance is set too low or 
                min_inter_distance is set too high, then less configurations may be generated.
        configurations_path - Local path to the file in to write the configurations to.
        min_distance        - The minimum distance between the centers of mass of the two molecules in any
                configuration.
        max_distance        - The maximum distance between the centers of mass of the two molecules in any
                configuration.
        min_inter_distance  - Minimum intermolecular distance in any config is the sum of two atoms vdw radii * this
                value.
        progression         - If True, use a linear progression for distance between atoms, otherwise random
                progression.r
        use_grid            - If False, configurations are space roughly evenly between min_distance and max_distance.
                If True, then configurations are placed at intervals along this distance based on step_size.
        step_size           - If use_grid is True, then this dictates the distance of the spacing interval used to
                place the centers of masses of the molecules. Otherwise, this parameter has no effect.
        num_attempts        - The number of tries to generate a config at any given distance before giving up and
                moving to the next one.
        logarithmic         - If True, then a logarithmic progression is used to generate the configurations.
                This means more configs are generated at lower distances.
        distribution        - An implementation of DistributionFunction. If specified, the logarithmic argument
                is ignored and this distribution is used to choose the distances between configurations. Should
                be implemented over the domain [0,1]. So the first config will have distance
                distribution.get_value(0) and the last config will have distance distribution.get_value(1).
        mol1_atom_index     - If specified, then the first molecule will be centered around the atom at this index
                rather than its center of mass.
        mol2_atom_index     - If specified, then the second molecule will be centered around the atom at this index
                rather than its center of mass.
        seed                - The same seed will generate the same configurations.

    Returns:
        None.
    """

    if mol1_atom_index is not None and mol2_atom_index is not None:
        config_generator = configurations.AtomDistanceConfigurationGenerator(settings_path, mol1_atom_index,
                                                                             mol2_atom_index,
                                                                             min_distance=min_distance,
                                                                             max_distance=max_distance,
                                                                             min_inter_distance=min_inter_distance,
                                                                             progression=progression,
                                                                             use_grid=use_grid, step_size=step_size,
                                                                             num_attempts=num_attempts,
                                                                             logarithmic=logarithmic,
                                                                             distribution=distribution)
    elif mol1_atom_index is None and mol2_atom_index is None:
        config_generator = configurations.DistanceSamplingConfigurationGenerator(settings_path,
                                                                                 min_distance=min_distance,
                                                                                 max_distance=max_distance,
                                                                                 min_inter_distance=min_inter_distance,
                                                                                 progression=progression,
                                                                                 use_grid=use_grid,
                                                                                 step_size=step_size,
                                                                                 num_attempts=num_attempts,
                                                                                 logarithmic=logarithmic,
                                                                                 distribution=distribution)
    else:
        raise InconsistentValueError("mol1_atom_index", "mol2_atom_index", mol1_atom_index, mol2_atom_index,
                                     "either one or both of these values must be specified")

    configurations.ConfigurationGenerator.generate_configs_from_file_to_file([geo1_path, geo2_path],
                                                                             configurations_path,
                                                                             config_generator,
                                                                             number_of_configs,
                                                                             seed=seed)
def generate_atom_distance_configurations(settings_path, geo1_path, geo2_path, number_of_configs, configurations_path,
        mol1_atom_index, mol2_atom_index, min_distance=1, max_distance=5, min_inter_distance=0.8, progression=False,
        use_grid=False, step_size=0.5, num_attempts=100, logarithmic=False, distribution=None, seed=None):
    """
    Generates 2b configurations for a given dimer by placing two atoms a certain distance apart and applying
    random rotations.

    Args:
        settings_path       - Local path to '.ini' settings file with all relevant settings.
        min_distance        - The minimum distance between the centers of mass of the two monomers.
        max_distance        - The maximum distance between the centers of mass of the two monomers.
        min_inter_distance  - Minimum intermolecular distance is this times the sum of the van der walls radii of two
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
        seed                - The same seed will generate the same configurations.

    Returns:
        None.
    """

    config_generator = configurations.AtomDistanceConfigurationGenerator(settings_path, mol1_atom_index, mol2_atom_index,
                                                                         min_distance=min_distance,
                                                                         max_distance=max_distance,
                                                                         min_inter_distance=min_inter_distance,
                                                                         progression=progression,
                                                                         use_grid=use_grid, step_size=step_size,
                                                                         num_attempts=num_attempts,
                                                                         logarithmic=logarithmic,
                                                                         distribution=distribution)

    configurations.ConfigurationGenerator.generate_configs_from_file_to_file([geo1_path, geo2_path],
                                                                             configurations_path,
                                                                             config_generator,
                                                                             number_of_configs,
                                                                             seed=seed)

def generate_configurations(settings_path, number_of_configs, configurations_path, *geo_paths, radius=10,
                            min_inter_distance=0.8, num_attempts=100, seed=None, logarithmic=False, distribution=None):
    """
    Generates a set of n body configurations by randomly placing monomer geometries in a sphere.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        number_of_configs   - Number of configurations to generate.
        configurations_path         - Local path to the file to write the configurations.
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
        distribution        - An implementation of DistributionFunction. If specified, the logarithmic argument
                is ignored and this distribution is used to choose the distances between configurations. Should
                be implemented over the domain [0,1]. So the first config will have distance =
                distribution.get_value(0) and the last config will have distance = distribution.get_value(1).

    Returns:
        None
    """

    config_generator = configurations.RandomSamplingConfigurationGenerator(settings_path,
                                                                           radius=radius,
                                                                           min_inter_distance=min_inter_distance,
                                                                           num_attempts=num_attempts,
                                                                           logarithmic=logarithmic,
                                                                           distribution=distribution)

    configurations.ConfigurationGenerator.generate_configs_from_file_to_file(geo_paths,
                                                                             configurations_path,
                                                                             config_generator,
                                                                             number_of_configs,
                                                                             seed=seed)

def init_database(settings_path, database_config_path, configurations_path, method, basis, cp, *tags, optimized = False):
    """
    Creates a database from the given configuration .xyz files. Can be called on a new database
    to create a new database, or an existing database to add more energies to be calculated

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        configurations_path - Local path to a single .xyz file.
        method              - QM method to use to calculate the energy of these configurations.
        basis               - QM basis to use to calculate the energy of these configurations.
        cp                  - Use counterpoise correction for these configurations?
        tags                - Mark the new configurations with these tags.
        optimized           - Are these configurations optimized geometries? Defualt is False.

    Returns:
        None.
    """

    database.initialize_database(settings_path, database_config_path, configurations_path, method, basis, cp, *tags, optimized = optimized)


def fill_database(settings_path, database_config_path, client_name, *tags, calculation_count = sys.maxsize, qm_options={}):
    """
    Goes through all the uncalculated energies in a database and calculates them. Will take a while. May be interrupted
    and restarted.
    
    Args:
        settings_path       - Local path to the file containing all relevant settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        client_name         - Name of the client performing these calculations.
        tags                - Only perform calculations marked with at least one of these tags.
        calculation_count   - Maximum number of calculations to perform. Unlimited if None.
        qm_options           - Dictionary of extra arguments to be passed to the QM code doing the calculation.

    Returns:
        None.
    """

    if calculation_count is None:
        calculation_count = sys.maxsize

    database.fill_database(settings_path, database_config_path, client_name, *tags, calculation_count=calculation_count, qm_options=qm_options)

def make_jobs(settings_path, database_config_path, client_name, job_dir, *tags, num_jobs=sys.maxsize, qm_options={}):
    """
    Makes a Job file for each energy that still needs to be calculated in this Database.

    Args:
        settings_path       - Local path to the file containing all relevant settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        client_name         - Name of the client that will perform these jobs
        job_dir             - Local path to the directory to place the job files in.
        tags                - Onlt  make jobs for calculations marked with at least one of these tags.
        num_jobs            - The number of jobs to generate. Unlimted if None.
        qm_options           - Dictionary of extra arguments to be passed to the QM code doing the calculation.

    Returns:
        None.
    """

    job_handler = database.get_job_handler(settings_path)

    job_handler.make_all_jobs(database_config_path, client_name, job_dir, *tags, num_jobs=num_jobs, qm_options=qm_options)

def read_jobs(settings_path, database_config_path, job_dir, overwrite=False):
    """
    Searches the given directory for completed job directories and enters
    the results into the database.

    Any directory that starts with job_ will be considered a completed job directory.
    After data is entered into the database, these directories will be renamed from
    job_* to job_*_done.

    Args:
        settings_path       - Local path to the file containing all relevant settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        job_dir             - Local path the the directory to search.

    Returns:
        None.
    """

    job_handler = database.get_job_handler(settings_path)

    job_handler.read_all_jobs(database_config_path, job_dir, overwrite=overwrite)

def generate_training_set(settings_path, database_config_path, training_set_path, method, basis,
        cp, *tags, e_bind_min=-float('inf'), e_bind_max=float('inf'), e_mon_min=-float('inf'), e_mon_max=float('inf'),
        deprecated_fitcode=False):
    """"
    Creates a training set file from the calculated energies in a database.

    Args:
        settings_path       - Local path to the ".ini" file with all relevent settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        training_set_path   - Local path to file to write training set to.
        method              - Use energies calculated with this method. Use % for any method.
        basis               - Use energies calculated with this basis. Use % for any basis.
        cp                  - Use energies calculated with this cp. Use 0 for False, 1 for True, or % for any cp.
        tags                - Use energies marked with at least one of these tags. Use % for any tag.
        e_bind_min          - Minimum binding energy allowed, inclusive.
        e_bind_max          - Maximum binding energy allowed, exclusive.
        e_mon_max           - Minimum monomer deformation energy allowed, inclusive.
        e_mon_max           - Maximum monomer deformation energy allowed, exclusive.
        deprecated_fitcode  - Is this function being called to be used with the deprecated fitcode?
                The output of the 1b and 2b training sets will be different.

    Return:
        None.
    """
    database.generate_training_set(settings_path, database_config_path, training_set_path, method, basis,
            cp, *tags, e_bind_min=e_bind_min, e_bind_max=e_bind_max, e_mon_min=e_mon_min, e_mon_max=e_mon_max,
            deprecated_fitcode=deprecated_fitcode)

def generate_poly_input(settings_path, molecule_in, in_file_path, virtual_sites=["X", "Y", "Z"]):
    """
    Generates an input file for polynomial generation.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        molecule_in         - String idicating symmetry of molecule, ie "A1B2_A1B2" for (CO2)2
        in_file_path        - Local path to the file to write the polynomial input to.
        virtual_sites       - List of Symmetry labels that are virtual sites.
                Default: ["X", "Y", "Z"]

    Returns:
        None.
    """

    polynomials.generate_input_poly(settings_path, molecule_in, in_file_path, virtual_sites=virtual_sites)

def generate_polynomials(settings_path, poly_in_path, order, poly_dir_path, generate_direct_gradients=True,
                         num_gradient_terms_per_line=1):
    """
    Generates polynomial input for maple and some ".cpp" and ".h" polynomial files.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        poly_in_path        - Local path to the file to read the polynomial input from. Name of file should be in
                the format "A1B2.in". It is ok to have extra directories prior to file name. For example:
                "thisplace/thatplace/A3.in".
        order               - The order of the polynomial to generate.
        poly_dir_path       - Local path to the directory to write the polynomial files in.
        generate_direct_gradients - If True, then a gradients cpp file is generate additionally to the polynomial
                cpp file. This may take a while.
                defualt: False.
        num_gradient_terms_per_line - The number of terms to put in each line in the cpp file
                for computing the gradients. Larger numbers may result in slower compilation times.

    Returns:
        None.
    """

    poly_generator = polynomials.PolynomialGenerator(settings_path,
                                                     num_gradient_terms_per_line=num_gradient_terms_per_line)

    poly_generator.generate_polynomial(poly_in_path, order, poly_dir_path,
            generate_direct_gradients=generate_direct_gradients)

def execute_maple(settings_path, poly_dir_path):
    """
    Runs maple on the maple input files generated during polynomial generation to convert them to C++ files.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        poly_directory      - Local path to the directory to read the ".maple" files from and write the ".cpp" files
                to.

    Returns:
        None.
    """

    system.format_print("Executing maple to create c polynomial files...", bold=True, color=system.Color.YELLOW)

    this_file_path = os.path.dirname(os.path.abspath(__file__))

    original_dir = os.getcwd()

    os.chdir(poly_dir_path)

    # clear any old files, because for some reason maple appends to existing files instead of clobbering
    with contextlib.suppress(FileNotFoundError):
        os.remove("poly-grd.c")
        os.remove("poly-nogrd.c")
        os.remove("poly-grd.cpp")
        os.remove("poly-nogrd.cpp")
    
    system.call("maple", "poly-grd.maple")
    system.call("maple", "poly-nogrd.maple")

    system.format_print("Maple executed successfully!", italics=True)
    system.format_print("Converting c polynomial files to cpp polynomial files", italics=True)

    clean_maple_exe = os.path.join(this_file_path, "polynomials", "scripts", "clean-maple-c.pl")

    with open("poly-grd.c", "r") as in_file, open("poly-grd.cpp", "w") as out_file:
        system.call(clean_maple_exe, in_file = in_file, out_file = out_file)
    with open("poly-nogrd.c", "r") as in_file, open("poly-nogrd.cpp", "w") as out_file:
        system.call(clean_maple_exe, in_file = in_file, out_file = out_file)

    system.format_print("cpp files generated successfully!", bold=True, color=system.Color.GREEN)

    os.chdir(original_dir)

def get_system_properties(settings_path, config_path, geo_paths,
                          distance_between=20,
                          use_published_polarizabilities=True,
                          method="wb97m-v",
                          basis="aug-cc-pvtz",
                          num_digits=4,
                          virtual_sites=["X", "Y", "Z"],use_cm5=True):
    """
        Obtains information such as charges and pols that will be needed for the fitting.

        Qchem is required for this step to work.

        Args:
            settings_path       - Local path to the file containing all relevent settings information.
            config_path         - Local path to file to write the config file to.
            geo_paths           - List of local paths to the optimized geometries to include in this fit config.
            distance_between    - The Distance between each geometry in the qchem calculation. If the qchem calculation
                    does not converge, try different values of this.
            use_published_polarizabilities - use published polarizabilites from
                    DOI: 10.1080/00268976.2018.1535143 rather than the ones Marc gave me to use.
                    Default: True.
            method              - Method to use for charges, polarizabilities, and c6 constants.
                    Default: wb97m-v.
            basis               - Basis to use for charges, polarizabilites, and c6 constants.
                    Default: aug-cc-pvtz.
            num_digits            - Number of digits after the decimal point to include in charges, c6, and polarizabilites.
                    Default: 4
            virtual_sites       - List of Symmetry labels that are virtual sites.
                    Default: ["X", "Y", "Z"]
            use_cm5             - Use CM5 charges rather than ChElPG

        Returns:
            charges             - List with the charges of each monomer [[mon1],[mon2],..]
            pols                - List with the pols of each monomer [[mon1],[mon2],..]
            C6                  - List with the C6 constants
        """

    chg, pol, c6 = fitting.get_system_properties(settings_path, config_path, geo_paths,
                                  distance_between=distance_between,
                                  use_published_polarizabilities=use_published_polarizabilities,
                                  method=method,
                                  basis=basis,
                                  num_digits=num_digits,
                                  virtual_sites=virtual_sites,use_cm5=use_cm5)
    return chg, pol, c6

def write_config_file(settings_file, config_path, charges,
                      pols, geo_paths, C6 = [0.0], polfacs = None,
                      d6 = None, A = None,
                      kmin = 0.0, kmax = 50.0, dmin = 0.0, dmax = 50.0, 
                      bmin = 0.0, bmax = 10.0,
                      kmin_init = 1.0, kmax_init = 4.0, dmin_init = 1.0, dmax_init = 4.0,
                      bmin_init = 1.0, bmax_init = 4.0,
                      r_in=6.0, r_out=8.0,
                      energy_range = 20, alpha = 0.0005,
                      virtual_sites_label = ['X','Y','Z'],
                      var_intra = "exp", var_inter = "exp", var_virtual_sites = "coul"):
    """
        Writes the config file.

        Qchem is required for this step to work.

        Args:
            settings_path       - Local path to the file containing all relevent settings information.
            config_path         - Local path to file to write the config file to.
            charges             - List with the charges of each monomer [[mon1],[mon2],..]
            pols                - List with the pols of each monomer [[mon1],[mon2],..]
            geo_paths           - List of local paths to the optimized geometries to include in this fit config.
            C6                  - List with the C6 constants
                    Default: [0.0]
            polfacs             - List with the polarizability factors of each monomer [[mon1],[mon2],..]
                    Default: None
            d6                  - List with the d6 constants
                    Default: None
            A                   - List with the buckingham A parameters
                    Default:  None
            kmin                - Minimum value of k allowed while fitting
                    Default: 0.0
            kmax                - Maximum value of k allowed while fitting
                    Default: 50.0
            dmin                - Minimum value of d allowed while fitting
                    Default: 0.0
            dmax                - Maximum value of d allowed while fitting
                    Default: 50.0
            bmin                - Minimum value of b allowed while fitting
                    Default: 0.0
            bmax                - Maximum value of b allowed while fitting
                    Default: 10.0
            kmin_init           - Minimum value of k allowed in initialization
                    Default: 1.0
            kmax_init           - Maximum value of k allowed in initialization
                    Default: 4.0
            dmin_init           - Minimum value of d allowed in initialization
                    Default: 1.0
            dmax_init           - Maximum value of d allowed in initialization
                    Default: 4.0
            bmin_init           - Minimum value of b allowed in initialization
                    Default: 1.0
            bmax_init           - Maximum value of b allowed in initialization
                    Default: 4.0
            r_in                - Distance at which polynomials start to decay to 0. Ideally, 
                                  it should be 2-3A beyond the peak of the first solvation shell.
                    Default: 6.0
            r_out               - Distance at which polynomials are 0. As a rule of thumb,
                                  it should be 2A beyond the inner cutoff.
                    Default: 8.0
            energy_range        - Value of DE in the weight expressions: w = (DE/(E-Emin+DE))^2
                    Default: 20.0
            alpha               - Ridge regression parameter
                    Default: 0.0005
            virtual_sites_label - List of Symmetry labels that are virtual sites.
                    Default: ["X", "Y", "Z"]
            var_intra           - Type of variable used for intramolecular distances.
                                  exp = exp(-kr)
                                  exp0 = exp(-k(r-r0))
                                  coul = exp(-kr)/r
                                  coul0 = exp(-k(r-r0))/r
                    Default: exp
            var_inter           - Type of variable used for intermolecular distances.
                                  exp = exp(-kr)
                                  exp0 = exp(-k(r-r0))
                                  coul = exp(-kr)/r
                                  coul0 = exp(-k(r-r0))/r
                    Default: exp
            var_virtual_sites   - Type of variable used for distances involving polynomial virtual sites.
                                  exp = exp(-kr)
                                  exp0 = exp(-k(r-r0))
                                  coul = exp(-kr)/r
                                  coul0 = exp(-k(r-r0))/r 
                    Default: coul

        Returns:
            None.
    """

    fitting.write_config_file(settings_file, config_path, charges,
                              pols, geo_paths, C6, polfacs,
                              d6, A,
                              kmin, kmax, dmin, dmax, bmin, bmax,
                              kmin_init, kmax_init, dmin_init, dmax_init,
                              bmin_init, bmax_init, r_in, r_out,
                              energy_range, alpha,
                              virtual_sites_label,
                              var_intra, var_inter, var_virtual_sites)




def generate_mbnrg_fitting_code(settings_path, config_path, poly_in_path, poly_path, poly_order, fit_dir_path, use_direct=False, version = "v1"):
    """
    Generates the fit code based on the polynomials for a system

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        config_path         - Local path to the dimer config file.
        poly_in_path        - Local path to the the A3B2.in type file to read polynomial input from.
        poly_path           - Local path to directory where polynomial files are.
        poly_order          - The order of the polynomial in poly_path.
        fit_dir_path        - Local path to directory to generate fit code in.
        use_direct          - If true, it will use the direct polynomials instead of the mapleoptimized

    Returns:
        None.
    """

    files.init_directory(fit_dir_path)

    if not os.path.isdir(fit_dir_path):
        os.mkdir(fit_dir_path)

    fitting.prepare_mbnrg_fitting_code(settings_path, config_path, poly_in_path, poly_path, poly_order, fit_dir_path, use_direct,version)

def generate_ttmnrg_fitting_code(settings_path, config_path, fit_dir_path):
    """
    Generates the TTM fit code based on settings and system configuration

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        config_path         - Local path to the dimer config file.
        fit_dir_path        - Local path to directory to generate fit code in.

    Returns:
        None.
    """

    files.init_directory(fit_dir_path)

    if not os.path.isdir(fit_dir_path):
        os.mkdir(fit_dir_path)

    fitting.prepare_ttmnrg_fitting_code(settings_path, config_path, fit_dir_path)

def compile_fit_code(settings_path, fit_dir_path):
    """
    Compiles the fit code in the given directory.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fit_dir_path        - Local path to the directory to read uncompiled fit code from and write compiled fit code
                to.

    Returns:
        None.
    """

    system.format_print("Compiling fit code...", bold=True, color=system.Color.YELLOW)

    original_dir = os.getcwd()

    os.chdir(fit_dir_path + "/src/")

    system.call("make", "clean")
    system.call("make")

    os.chdir(original_dir)

    system.format_print("Fit code compilation successful!", bold=True, color=system.Color.GREEN)

#def add_A_and_b_to_config_file(settings_path, fit_dir_path, config_path):
#    
#    # read d6 and A constants from ttm output file
#    with open(os.path.join(fit_dir_path, "ttm-params.txt"), "r") as ttm_file:
#        A = [float(a) for a in ttm_file.readline().split()]
#        d6 = [float(d) for d in ttm_file.readline().split()]
#
#    # write d6 and A to the config.ini file
#    lines = []
#    with open(config_path, "r") as config_file:
#        for line in config_file:
#            lines.append(line)
#
#    found_A = False
#    with open(config_path, "w") as config_file:
#        for line in lines:
#            if line.startswith("A = "):
#                found_A = True
#            if not line.startswith("d6 = "):
#                config_file.write(line)
#
#            else:
#                config_file.write("d6 = {}\n".format([[], [], d6]))
#
#        if not found_A:
#            config_file.write("A = {}\n".format([[], [], A]))


def prepare_fits(settings_path, fit_dir_path, training_set_path, fits_path, DE=20, alpha=0.0005, num_fits=10, ttm=False,
                 over_ttm=False):
    """
    Prepares fits to be run by creating directories with executable scripts in them.

    Each directory will contain a run_fit.sh script, which when ran will run one fit.

    The directories will appear in the directory with your other log files.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fit_dir_path        - Local path to the directory containing the compiled fitcode.
        training_set_path   - Local path to training set to use for all the fits.
        fits_path           - Local path to the directory to create the fits in.
        DE                  - Low DE places more weight on low energy training set items. 
                              Large DE places even weight on all training set items. Weights w_n are computed as
                              w_n = (DE / (E_n - E_min + DE))**2
        alpha               - Weight for the regularization parameter in the fits. Large alpha means coefficients will
                be smaller. (Less likely for one or more coefficients to blow up.)
        num_fits            - How many new directories with executables to make. Existing ones will not be changed.
        ttm                 - True if these are ttm fits. False otherwise.
        over_ttm            - Only used if ttm is False, if enabled, will fit polynomials over ttm.

    Returns:
        None.
    """
    # Get information
    settings = SettingsReader(settings_path)
    nb = len(settings.get("molecule","names").split(","))
    workdir = os.getcwd()

    fit_folder_prefix = os.path.join(workdir, fits_path)
    if not os.path.exists(fit_folder_prefix):
        os.mkdir(fit_folder_prefix)

    if ttm:
        fit_executable_path = os.path.join(workdir, fit_dir_path, "bin/fit-{}b-ttm".format(nb))
    elif over_ttm:
        fit_executable_path = os.path.join(workdir, fit_dir_path, "bin/fit-{}b-over-ttm".format(nb))
    else:
        fit_executable_path = os.path.join(workdir, fit_dir_path, "bin/fit-{}b".format(nb))

    # initialize indexes
    fit_index = 1
    count = 0
    # Prepare fits
    while count < num_fits:
        os.chdir(fit_folder_prefix)
        new_fit_folder = "fit" + str(fit_index)
        if os.path.exists(new_fit_folder):
            print("{} folder already exists. Trying next index.".format(new_fit_folder))
            fit_index += 1
        else:
            os.mkdir(new_fit_folder)
            os.chdir(new_fit_folder)

            training_set_basename = os.path.basename(training_set_path)

            # Link the training set to the fit folder
            if training_set_path.startswith('/'):
                # we have an absolute path
                system.call("ln", "-s", training_set_path, training_set_basename)
            else:
                # this is a relative path to working directory
                system.call("ln", "-s", "{}/{}".format(workdir, training_set_path), training_set_basename)

            # Create bash script that will run the fit
            with open("run_fit.sh",'w') as my_bash:
                my_bash.write("""#!/bin/bash
restarted=0
if [ -f "fit-{nb}b-restart-1.cdl" ]
then
    ncgen -o restart.nc fit-{nb}b-restart-1.cdl
    if [ $? -ne 0 ]
    then
        echo "Cannot restart from fit-{nb}b-restart-1.cdl, it seems to be corrupted."
    else
        restarted=1
        echo "Restarting fit from fit-{nb}b-restart-1.cdl."
    fi
fi

if [ -f "fit-{nb}b-restart-2.cdl" ] && [ $restarted -eq 0 ]
then
    ncgen -o restart.nc fit-{nb}b-restart-2.cdl
    if [ $? -ne 0 ]
    then
        echo "Cannot restart from fit-{nb}b-restart-2.cdl, it seems to be corrupted."
    else
        restarted=1
        echo "Restarting fit from fit-{nb}b-restart-2.cdl."
    fi
fi


if [ $restarted -eq 1 ]
then
    {fit_exe} {train_set} {DE} {alpha} restart.nc > fit.log 2> fit.err 
else
    {fit_exe} {train_set} {DE} {alpha} > fit.log 2> fit.err 
fi
""".format(**{"nb": nb, "fit_exe": fit_executable_path, "train_set": training_set_basename, "DE": DE, "alpha": alpha}))

            system.call("chmod", "744", "run_fit.sh")
            fit_index += 1
            count += 1 
            print("Succesfully created fit folder {}.".format(new_fit_folder))

    os.chdir(workdir)

def execute_fits(settings_path, fits_path):
    """
    Looks for fit executables in fit directories in the logs directory and runs the fits.

    Will ignore any fits that already have a log file - suggesting they have already been run.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fits_path           - Local path to the directory to create the fits in.

    Returns:
        None.
    """

    settings = SettingsReader(settings_path)
    nb = len(settings.get("molecule","names").split(","))

    # Get information
    workdir = os.getcwd()
    fit_folder_prefix = os.path.join(workdir, fits_path)

    # A fit folder that does not have a fit.log file inside is considered not run
    # In that case, the run_fit.sh will be executed
    os.chdir(fit_folder_prefix)
    all_fits = glob.glob("fit*")
    for fit in all_fits:
        os.chdir(fit)
        if not os.path.exists("fit-{}b.cdl".format(nb)):
            print("{} is running.".format(fit))
            print(os.path.dirname(os.curdir))
            system.call("./run_fit.sh", out_file = None)
            print("{} is completed.".format(fit))
        else:
            print("{} is already done. Continuing...".format(fit))
        os.chdir("../")

    os.chdir(workdir)
    
def retrieve_best_fit(settings_path, fits_path, fitted_nc_path = "mbnrg.nc", fitted_ttmnrg_params = "ttm-nrg_params.dat"):
    """
    Looks through all log files in all fit directories in the log path and finds the best fit.

    The best_fit will end up inside a directory in the logs directory called "best_fit".

    If any new fit is better than the current best fit, it will replace best_fit with the new
    best one.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fits_path           - Local path to the directory to create the fits in.
        fitted_nc_path      - Generate a .nc file with the parameters for the best fit at this location.
        fitted_ttmnrg_params - Rename the output where the TTM-nrg params are to this name

    Returns:
        None.
    """
    # Get information
    settings = SettingsReader(settings_path)
    workdir = os.getcwd()
    fit_folder_prefix = os.path.join(workdir, fits_path)

    # Loop over all the fits, check the output, and store the results.
    os.chdir(fit_folder_prefix)
    all_fits = glob.glob("fit*")
    results = []
    for fit in all_fits:
        os.chdir(fit)

        try:
            with open("fit.log",'r') as logfile:
                log_lines = logfile.readlines()
                converged = "!!! converged !!!\n" in log_lines
                full_rmsd = float(log_lines[-7].split()[2])
                wfull_rmsd = float(log_lines[-6].split()[2])
                max_error = float(log_lines[-5].split()[2])
                low_rmsd = float(log_lines[-4].split()[2])
                low_max_error = float(log_lines[-3].split()[2])

                if not converged:
                    print("Looks like the fit in " + fit + " didn't converge...")
                    print("Maybe you want to rerun " + fit + " again.")
                    results.append([fit, full_rmsd, wfull_rmsd, max_error, low_rmsd, low_max_error, False])
                else:
                    results.append([fit, full_rmsd, wfull_rmsd, max_error, low_rmsd, low_max_error, True])
        except (FileNotFoundError, IndexError, ValueError):
            print("Doesn't seem that the log file in " + fit + " is correct...")
            print("Maybe you want to rerun " + fit + " again.")
            results.append([fit, float('inf'), float('inf'), float('inf'), float('inf'), float('inf'), False])

        os.chdir("../")

    # Sort the results according to weighet RMSD of the full TS
    sorted_results = sorted(results, key=lambda x: x[2])

    # Get the best result
    best_fit = sorted_results[0][0]
    best_results = sorted_results[0]

    print("Best fit is fit {}.".format(best_fit))

    if not best_results[6]:
        print("Fit with lowest RMSD did not converge, maybe run it again with more iterations?")
        print("Still using that fit as the best one.")
    
    # Check if best fit folder exists.
    # If it does exist, check the output in best folder to ensure that the new best
    # fit is actually the best
    # If not, just store the best fit in there
    if not os.path.exists("best_fit"):
        print("No previous best fit, so just storing the current best fit.")
        system.call("ln", "-s", best_fit, "best_fit")
    else:
        os.chdir("best_fit")
        with open("fit.log",'r') as logfile:
            log_lines = logfile.readlines()
            converged = "!!! converged !!!\n" in log_lines
            full_rmsd = float(log_lines[-7].split()[2])
            wfull_rmsd = float(log_lines[-6].split()[2])
            max_error = float(log_lines[-5].split()[2])
            low_rmsd = float(log_lines[-4].split()[2])
            low_max_error = float(log_lines[-3].split()[2])
            previous_best = [fit, full_rmsd, wfull_rmsd, max_error, low_rmsd, low_max_error, converged]
        os.chdir("../")

        replace = previous_best[2] > best_results[2]

        if replace:
            print("New best fit is better than the old one.")
            print("Replacing best_fit by {}.".format(sorted_results[0][0]))
            system.call("rm", "-rf", "best_fit")
            system.call("cp", "-r", sorted_results[0][0], "best_fit")
        else:
            print("New best fit is not better than the old one.")
            print("Keeping the old best fit.".format(sorted_results[0][0]))
            best_results = previous_best


    os.chdir("best_fit")
    nb = len(settings.get("molecule","names").split(","))
    if os.path.exists("fit-" + str(nb) + "b.cdl"):
        system.call("ncgen", "-o", fitted_nc_path, "fit-" + str(nb) + "b.cdl")
    elif os.path.exists("ttm-nrg_params.dat") and fitted_ttmnrg_params != "ttm-nrg_params.dat":
        os.system("mv ttm-nrg_params.dat " + fitted_ttmnrg_params)

    # Report best RMSD
    print("Best fit found has a weighted RMSD of {} kcal/mol, a low energy RMSD of {} kcal/mol, and a maximum error in the low energy training set of {} kcal/mol".format(best_results[2], best_results[4], best_results[5]))

    os.chdir(workdir)

def update_config_with_ttm(settings_path, fits_path, config_path, fitted_ttmnrg_params = "ttm-nrg_params.dat"):
    """
    Updates a fit config.ini file with the A and d constants from the best ttm fit.

    The best ttm fit values will be taken from the best_fit directory inside the logs directory.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fits_path           - Local path to the directory to create the fits in.
        config_path         - Local path to the config file to update.
        fitted_ttmnrg_params - Name of the output where the TTM-nrg params


    Returns:
        None.
    """
    config = SettingsReader(config_path)
    workdir = os.getcwd()
    fit_folder_prefix = os.path.join(workdir, fits_path)

    with open(fit_folder_prefix + "/best_fit/" + fitted_ttmnrg_params, 'r') as ttm_file:
        a_buck = [float(a) for a in ttm_file.readline().strip().split()]
        b_buck = [float(b) for b in ttm_file.readline().strip().split()]

    config.set("fitting", "A", a_buck)
    config.set("fitting", "d6", b_buck)

    config.write(config_path)


def generate_MBX_files(settings_path, config_file, mon_ids, degree, do_ttmnrg = False, mbnrg_fits_path = None, MBX_HOME=None, version="v1",
                            virtual_sites=["X", "Y", "Z"]):
    """
    Creates the pieces of code that will be needed in MBX to use this potential

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        config_file         - Local path to the config file to update.
        mon_ids             - List with the monomer ids that will be used in MBX to identify the monomers added
        do_ttmnrg           - If true, TTM-nrg will also be added. Should be used only when doing a 2b fit.
        mbnrg_fits_path     - Path to the folder that contains the mb-nrg fits. It is expected to have a "best_fit" folder inside that one.
        degree              - Order of the polynomial (if any)
        MBX_HOME            - If defined, the pieces of code will be automatically added to MBX, and the necessary files will also be copied
        version             - iString with the extension to the polynomial and namespaces that will be added to differentiate between different versions of the same monomer fit.
        virtual_sites       - List of the caracters that represent virtual sites (not real atoms) used in the polynomials (i.e., lone pairs)

    Returns:
        None.
    """

    fitting.generate_software_files(settings_path, config_file, mon_ids, do_ttmnrg, mbnrg_fits_path, degree, MBX_HOME, version, virtual_sites)

def calculate_model_energies(settings_path, fitting_code_dir_path, fits_path, configurations,
                             ttm=False, over_ttm=False, nc_path = "mbnrg.nc",
                             fitted_ttmnrg_params = "ttm-nrg_params.dat"):
    """
    Returns a list with the energies of each configuration in configurations

    Args:
        settings_path         - Local path to the file containing all relevent settings information.
        fitting_code_dir_path - Local path to the directory containing the compiled fitcode.
        fits_path             - Local path to folder containing the fits
        configurations        - Configurations to evaluate.
        ttm                   - True if these are ttm fits. False otherwise.
        over_ttm              - Only used if ttm is False, if enabled, will fit polynomials over ttm.
        nc_path               - Netcdf file with the parameters for the best fit.
        fitted_ttmnrg_params  - Name of the output where the TTM-nrg params
    """

    # Get information
    workdir = os.getcwd()
    settings = SettingsReader(settings_path)
    nb = len(settings.get("molecule","names").split(","))

    if ttm:
        path_to_eval = os.path.join(workdir, fitting_code_dir_path, "bin/eval-{}b-ttm".format(nb))
        path_to_params = os.path.join(workdir, fits_path, "best_fit",fitted_ttmnrg_params)
    elif over_ttm:
        path_to_eval = os.path.join(workdir, fitting_code_dir_path, "bin/eval-{}b-over-ttm".format(nb))
        path_to_params = os.path.join(workdir, fits_path, "best_fit", nc_path)
    else:
        path_to_eval = os.path.join(workdir, fitting_code_dir_path, "bin/eval-{}b".format(nb))
        path_to_params = os.path.join(workdir, fits_path, "best_fit", nc_path)


    eval_obj = fitting.Evaluator(settings, path_to_eval)

    fit_energies = eval_obj.calculate_energies(path_to_params, configurations, is_training_format = False)
    return fit_energies

def get_correlation_data(settings_path, fitting_code_dir_path, fits_path, training_set,
                         split_energy = None, 
                         min_energy_plot = 0.0, max_energy_plot = 50.0,
                         correlation_prefix = "correlation",
                         correlation_directory = "correlation", minor_tick = 5.0, colors = None,
                         labels = None,
                         ttm=False, over_ttm=False, nc_path = "mbnrg.nc",
                         fitted_ttmnrg_params = "ttm-nrg_params.dat", max_1b = 100.0):
    """
    Generates correltation data for the training/test set passed as argument

    Args:
        settings_path         - Local path to the file containing all relevent settings information.
        fitting_code_dir_path - Local path to the directory containing the compiled fitcode.
        fits_path             - Local path to folder containing the fits
        training_set          - Configurations to evaluate. They need the binding energy and n-b energy in this order in the comment line.
        split_energy          - If not None, will split the correlation plot in two sets, low and high, depending on the value of the variable.
        min_energy_plot       - Lower bound of the energy in the plot
        max_energy_plot       - Upper bound of the energy in the plot
        correlation_prefix    - Prefix for the correlation files that will be generated.
        correlation_directory - Directory where all the correlation files will be put.
        minor_tick            - Interval of the minor ticks in the plot
        colors                - List of two elements with the colors for low energy and high energy in the correlation plot
        labels                - List of two elements with the x-axis and y-axis labels
        ttm                   - True if these are ttm fits. False otherwise.
        over_ttm              - Only used if ttm is False, if enabled, will fit polynomials over ttm.
        nc_path               - Netcdf file with the parameters for the best fit.
        fitted_ttmnrg_params  - Name of the output where the TTM-nrg params
        max_1b                - Maximum deformation energy allowed in the plot (IE - BE)
    """

    # Get information
    workdir = os.getcwd()
    settings = SettingsReader(settings_path)
    nb = len(settings.get("molecule","names").split(","))

    corr_folder_prefix = os.path.join(workdir, correlation_directory)
    if not os.path.exists(corr_folder_prefix):
        os.mkdir(corr_folder_prefix)

    if ttm:
        path_to_eval = os.path.join(workdir, fitting_code_dir_path, "bin/eval-{}b-ttm".format(nb))
        path_to_params = os.path.join(workdir, fits_path, "best_fit",fitted_ttmnrg_params)
    elif over_ttm:
        path_to_eval = os.path.join(workdir, fitting_code_dir_path, "bin/eval-{}b-over-ttm".format(nb))
        path_to_params = os.path.join(workdir, fits_path, "best_fit", nc_path)
    else:
        path_to_eval = os.path.join(workdir, fitting_code_dir_path, "bin/eval-{}b".format(nb))
        path_to_params = os.path.join(workdir, fits_path, "best_fit", nc_path)


    eval_obj = fitting.Evaluator(settings, path_to_eval)

    eval_obj.calculate_energies(path_to_params, training_set)

    correlation_file = correlation_prefix + ".dat"
    energies = eval_obj.write_correlation_file(correlation_file = correlation_file , split_energy = split_energy, max_1b = max_1b)

    os.system("mv *" +  correlation_file + " " + corr_folder_prefix)

    eval_obj.plot(do_ttm = ttm, split_energy = split_energy, correlation_prefix = correlation_prefix, min_e = min_energy_plot, max_e = max_energy_plot, minor_tick = minor_tick, colors = colors,labels = labels)

    os.system("mv *.png *.pdf " + " " + corr_folder_prefix)

    return energies


