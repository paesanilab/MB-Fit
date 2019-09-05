# external package imports
import os, sys, contextlib, glob

# local module imports
from .utils import SettingsReader, files, system, constants
from . import configurations, database, polynomials, fitting
from .database import Database
from .molecule import xyz_to_molecules


def apply_standard_order(settings_path, geo_path):

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


def optimize_geometry(settings_path, unopt_geo_path, opt_geo_path, method, basis):
    """
    Optimizes the geometry of the given molecule.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        unopt_geo_path      - Local path to the file to read the unoptimized geoemtry from.
        opt_geo_path        - Local path to the file to write the optimized geometry to.
        method              - The method to use for this geometry optimization.
        basis               - The basis to use for this geometry optimization.

    Returns:
        None.
    """

    configurations.optimize_geometry(settings_path, unopt_geo_path, opt_geo_path, method, basis)

def generate_normal_modes(settings_path, opt_geo_path, normal_modes_path, method, basis):
    """
    Generates the normal modes for the given molecule.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        opt_geo_path        - Local path to the file to read the optimized geometry from.
        normal_modes_path   - Local path to the file to write the normal modes to.
        method              - The method to use for this normal modes calculation.
        basis               - The basis to use for this normal modes calculation.

    Returns:
        Null dimension of normal modes.
    """
    
    dim_null = configurations.generate_normal_modes(settings_path, opt_geo_path, normal_modes_path, method, basis)

    return dim_null

def generate_normal_mode_configurations(settings_path, opt_geo_path, normal_modes_path, configurations_path,
        number_of_configs = 100, seed = None, temperature = None):
    """
    Generates normal mode configurations for a given monomer (or dimer or trimer) from a set of normal modes.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        opt_geo_path        - Local path to the file to read optimized geometry from.
        normal_modes_path   - Local path to the file to read normal modes from.
        dim_null            - The null dimension of this molecule, see generate_normal_modes().
        config_path         - Local path to the file to write configurations to.
        number_of_configs   - Number of configurations to generate
        seed                - The same seed with the same molecule and normal modes will always generate the same
                configurations.
        temperature         - Temperature at which normal mode sampling is done. If specified, configurations
                will use clasical normal mode distribution at the specified temperature instead of either geometric
                or linear progression.

    Returns:
        None.
    """

    if temperature is not None:
        temperature *= constants.kelvin_to_au
    configurations.generate_normal_mode_configurations(settings_path, opt_geo_path, normal_modes_path, configurations_path,
            number_of_configs, seed = seed, temperature=temperature)

def generate_2b_configurations(settings_path, geo1_path, geo2_path, number_of_configs, configurations_path, 
        min_distance = 1, max_distance = 5, min_inter_distance = 0.8, progression = False, use_grid = False, 
        step_size = 0.5, num_attempts = 100, logarithmic = False, seed = None):
    """
    Generates 2b configurations for a given dimer.

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
        seed                - The same seed will generate the same configurations.

    Returns:
        None.
    """

    configurations.generate_2b_configurations(settings_path, geo1_path, geo2_path, number_of_configs, configurations_path,
            min_distance, max_distance, min_inter_distance, progression, use_grid, step_size, num_attempts, logarithmic, seed)

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

    configurations.generate_configurations(settings_path, number_of_configs, config_path, *geo_paths, radius=radius, min_inter_distance=min_inter_distance, num_attempts=num_attempts,
                                      seed=seed, logarithmic=logarithmic)

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


def fill_database(settings_path, database_config_path, client_name, *tags, calculation_count = sys.maxsize):
    """
    Goes through all the uncalculated energies in a database and calculates them. Will take a while. May be interrupted
    and restarted.
    
    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        client_name         - Name of the client performing these calculations.
        tags                - Only perform calculations marked with at least one of these tags.
        calculation_count   - Maximum number of calculations to perform. Unlimited if None.

    Returns:
        None.
    """

    if calculation_count is None:
        calculation_count = sys.maxsize

    database.fill_database(settings_path, database_config_path, client_name, *tags, calculation_count=calculation_count)

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

def generate_1b_training_set(settings_path, database_config_path, training_set_path, molecule_name, method, basis, cp, *tags, e_min = 0, e_max = float('inf')):
    """
    Generates a 1b training set from the energies inside a database.

    ***deprecated, please use generate_training_set instead***

    Specific method, basis, and cp may be specified to only use energies calculated
    with a specific model.

    '%' can be used to stand in as a wild card, meaning any method/basis/cp will be used in the training set.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        training_set_path   - Local path to the file to write the training set to.
        molecule_name       - The name of the moelcule to generate a training set for.
        method              - Only use energies calcualted by this method.
        basis               - Only use energies calculated in this basis.
        cp                  - Only use energies calculated with this cp. Note that counterpoise correct has no
                effect on 1b energies.
        tags                - Only use energies marked with one or more of these tags.
        e_min               - Minimum (inclusive) energy of any config to include in this training set.
        e_max               - Maximum (exclusive) energy of any config to include in this training set.

    Returns:
        None.
    """

    database.generate_training_set(settings_path, database_config_path, training_set_path, method, basis,
        cp, *tags, e_bind_min=e_min, e_bind_max=e_max, e_mon_min=e_min, e_mon_max=e_max)


def generate_2b_training_set(settings_path, database_config_path, training_set_path, molecule_name, method, basis, cp, *tags,
            e_bind_max = float('inf'), e_mon_max = float('inf')):
    """
    Generates a 2b training set from the energies inside a database.

    ***deprecated, please use generate_training_set instead***

    Specific method, basis, and cp may be specified to only use energies calculated
    with a specific model.

    '%' can be used to stand in as a wild card, meaning any method/basis/cp will be used in the training set.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        training_set_path   - Local path to the file to write the training set to.
        molecule_name       - The name of the dimer.
        method              - Only use energies calculated by this method.
        basis               - Only use energies calculated in this basis.
        cp                  - Only use energies calculated with this cp.
        tags                - Only use energies marked with one or more of these tags.
        e_bind_max          - Maximum binding energy allowed
        e_mon_max           - Maximum monomer deformation energy allowed

    Returns:
        None.
    """
    database.generate_training_set(settings_path, database_config_path, training_set_path, method, basis,
        cp, *tags, e_bind_min=-float('inf'), e_bind_max=e_bind_max, e_mon_min=-float('inf'), e_mon_max=e_mon_max)

def generate_poly_input(settings_path, molecule_in, in_file_path):
    """
    Generates an input file for polynomial generation.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        molecule_in         - String idicating symmetry of molecule, ie "A1B2_A1B2" for (CO2)2
        in_file_path        - Local path to the file to write the polynomial input to.

    Returns:
        None.
    """

    polynomials.generate_input_poly(settings_path, molecule_in, in_file_path)

def generate_polynomials(settings_path, poly_in_path, order, poly_dir_path):
    """
    Generates polynomial input for maple and some ".cpp" and ".h" polynomial files.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        poly_in_path        - Local path to the file to read the polynomial input from. Name of file should be in
                the format "A1B2.in". It is ok to have extra directories prior to file name. For example:
                "thisplace/thatplace/A3.in".
        order               - The order of the polynomial to generate.
        poly_dir_path       - Local path to the directory to write the polynomial files in.

    Returns:
        None.
    """

    polynomials.generate_poly(settings_path, poly_in_path, order, poly_dir_path)

def execute_maple(settings_path, poly_dir_path):
    """
    Runs maple on the polynomial files in the specified directory to turn them into actual cpp files.

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

    with open("poly-grd.c", "r") as in_file, open("poly-grd.cpp", "w") as out_file:
        system.call("clean-maple-c.pl", in_file = in_file, out_file = out_file)
    with open("poly-nogrd.c", "r") as in_file, open("poly-nogrd.cpp", "w") as out_file:
        system.call("clean-maple-c.pl", in_file = in_file, out_file = out_file)

    system.format_print("cpp files generated successfully!", bold=True, color=system.Color.GREEN)

    os.chdir(original_dir)

def generate_fitting_config_file_new(settings_file, config_path, geo_paths, distance_between = 20, use_published_polarizabilities = True):
    """
        Generates the config file needed to perform a fit.

        Qchem is required for this step to work for 1 and 2 b.

        For 1B, a qchem calcualtion is performed and charges, polarizabilities, and c6 constants are read from the output.

        For 2B, a chem calculation is performed and intermolecular c6 cosntants are read from it.
        Charges, polarizabilities, and intramolecular c6 are read from the config_1b_paths.

        For 3B and above, charges, polarizabilities, and intramolecular c6 constants are read from the config_1b_paths.
        Intermolecular c6 constants are read from config_2b_paths.

        Args:
            settings_path       - Local path to the file containing all relevent settings information.
            config_path         - Local path to file to write the config file to.
            geo_paths           - List of local paths to the optimized geometries to include in this fit config.
            config_1b_paths     - List of local paths to 1b config files. Only used for 2B and above. Should be one
                    config for each monomer.
            config_2b_paths     - List of local paths to 2b config files. Only used for 3B and above. Should be one
                    config for each combination of monomers.
            distance_between    - The Distance between each geometry in the qchem calculation. If the qchem calculation
                    does not converge, try different values of this.
            use_published_polarizabilities - use published polarizabilites from
                    DOI: 10.1080/00268976.2018.1535143 rather than the ones Marc gave me to use.

        Returns:
            None.
        """
    fitting.generate_fitting_config_file_new(settings_file, config_path, geo_paths, distance_between=distance_between, use_published_polarizabilities=use_published_polarizabilities)


def generate_fitting_config_file(settings_file, config_path, geo_paths, config_1b_paths = [], config_2b_paths = [], distance_between = 20, use_published_polarizabilities = True):
    """
        Generates the config file needed to perform a fit.

        Qchem is required for this step to work for 1 and 2 b.

        For 1B, a qchem calcualtion is performed and charges, polarizabilities, and c6 constants are read from the output.

        For 2B, a chem calculation is performed and intermolecular c6 cosntants are read from it.
        Charges, polarizabilities, and intramolecular c6 are read from the config_1b_paths.

        For 3B and above, charges, polarizabilities, and intramolecular c6 constants are read from the config_1b_paths.
        Intermolecular c6 constants are read from config_2b_paths.

        Args:
            settings_path       - Local path to the file containing all relevent settings information.
            config_path         - Local path to file to write the config file to.
            geo_paths           - List of local paths to the optimized geometries to include in this fit config.
            config_1b_paths     - List of local paths to 1b config files. Only used for 2B and above. Should be one
                    config for each monomer.
            config_2b_paths     - List of local paths to 2b config files. Only used for 3B and above. Should be one
                    config for each combination of monomers.
            distance_between    - The Distance between each geometry in the qchem calculation. If the qchem calculation
                    does not converge, try different values of this.
            use_published_polarizabilities - use published polarizabilites from
                    DOI: 10.1080/00268976.2018.1535143 rather than the ones Marc gave me to use.

        Returns:
            None.
        """

    fitting.generate_fitting_config_file(settings_file, config_path, geo_paths, config_1b_paths=config_1b_paths,
            config_2b_paths=config_2b_paths, distance_between=distance_between, use_published_polarizabilities=use_published_polarizabilities)
    
def generate_1b_fit_code(settings_path, config_path, molecule_in, poly_in_path, poly_dir_path, order, fit_dir_path):
    """
    Generates the fit code based on the polynomials and config file for a monomer.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        config_path         - Local path to the monomer config file to read config info from.
        poly_in_path        - Local path to the file to read the polynomial input from. Name of file should be in
                the format A1B2.in, it is ok to have extra directories prior to file name (thisplace/thatplace/A3.in).
        poly_dir_path       - Local path to the directory where the polynomial ".h" and ".cpp" files are files are.
        order               - The order of the polynomial to in poly_dir_path.
        fit_dir_path        - Local path to the directory to write the fit code in.

    Returns:
        None.
    """

    fitting.prepare_1b_fitting_code(config_path, molecule_in, poly_in_path, poly_dir_path, order, fit_dir_path)

def generate_2b_ttm_fit_code(settings_path, config_path, molecule_in, fit_dir_path):
    """
    Generates the fit TTM fit code for a dimer.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        config_path         - Local path to the dimer config file to read config info from.
        molecule_in         - A String of fromat "A1B2". Same as poly_in_path but without ".in".
        fit_dir_path        - Local path to the directory to write the ttm fit code in.

    Returns:
        None.
    """

    this_file_path = os.path.dirname(os.path.abspath(__file__))

    original_dir = os.getcwd()

    files.init_directory(fit_dir_path)

    os.chdir(fit_dir_path)

    codes_path = os.path.join(this_file_path, "..", "codes", "2b-codes")
    
    # FOR SOME REASON THE SECOND LINE WORKS BUT NOT THE FIRST, WHY? WHO KNOWS.
    # system.call("cp", os.path.join(codes_path, "template/*"), ".")
    os.system("cp " + os.path.join(codes_path, "template", "*") + " .")

    ttm_script_path = os.path.join(codes_path, "get_2b_TTM_codes.py")

    if not settings_path.startswith("/"):
        settings_path = "{}/{}".format(original_dir, settings_path)
    if not config_path.startswith("/"):
        config_path = "{}/{}".format(original_dir, config_path)

    system.call("python", ttm_script_path, settings_path, config_path, molecule_in)

    os.chdir(original_dir)   
 
def generate_mbnrg_fitting_code(settings_path, config_path, poly_in_path, poly_path, poly_order, fit_dir_path):
    """
    Generates the fit code based on the polynomials for a system

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        config_path         - Local path to the dimer config file.
        poly_in_path        - Local path to the the A3B2.in type file to read polynomial input from.
        poly_path           - Local path to directory where polynomial files are.
        poly_order          - The order of the polynomial in poly_path.
        fit_dir_path        - Local path to directory to generate fit code in.

    Returns:
        None.
    """

    files.init_directory(fit_dir_path)

    if not os.path.isdir(fit_dir_path):
        os.mkdir(fit_dir_path)

    fitting.prepare_fitting_code(settings_path, config_path, poly_in_path, poly_path, poly_order, fit_dir_path)


def generate_2b_fit_code(settings_path, config_path, poly_in_path, poly_path, poly_order, fit_dir_path):
    """
    Generates the fit code based on the polynomials for a monomer

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        config_path         - Local path to the dimer config file.
        poly_in_path        - Local path to the the A3B2.in type file to read polynomial input from.
        poly_path           - Local path to directory where polynomial files are.
        poly_order          - The order of the polynomial in poly_path.
        fit_dir_path        - Local path to directory to generate fit code in.

    Returns:
        None.
    """

    files.init_directory(fit_dir_path)

    if not os.path.isdir(fit_dir_path):
        os.mkdir(fit_dir_path)

    fitting.prepare_2b_fitting_code(settings_path, config_path, poly_in_path, poly_path, poly_order, fit_dir_path)
 
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

    os.chdir(fit_dir_path)

    system.call("make", "clean")
    system.call("make")

    os.chdir(original_dir)

    system.format_print("Fit code compilation successful!", bold=True, color=system.Color.GREEN)

def perform_1b_fits(settings_path, fit_code_path, training_set_path, fit_dir_path, num_fits = 10):

    system.format_print("Performing {} fits from which the best will be chosen...".format(num_fits), bold=True, color=system.Color.YELLOW)

    settings = SettingsReader(settings_path)

    # init the required files
    files.init_directory(fit_dir_path)
    best_fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "mb", "best-fit.log"))
    fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "mb" "fit.log"))

    # keeps tracks of the number of attempts to generate a fit
    attempts = 1

    # generate an initial fit
    with open(best_fit_log_path, "w") as best_fit_log:
        system.call(fit_code_path, training_set_path, out_file = best_fit_log)
        os.rename("fit-1b.cdl", "best-fit-1b.cdl")
        os.rename("fit-1b-initial.cdl", "best-fit-1b-initial.cdl")
        os.rename("correlation.dat", "best-correlation.dat")

    with open(best_fit_log_path, "r") as best_fit_log:
        best_log_lines = best_fit_log.readlines()

    best_rmsd = float(best_log_lines[-6].split()[2])

    system.format_print("Completed first fit with rmsd {}.\n".format(best_rmsd), italics=True)

    while(attempts < num_fits):

        # generate a new fit
        with open(fit_log_path, "w") as fit_log:
            system.call(fit_code_path, training_set_path, out_file = fit_log)
  
        with open(fit_log_path, "r") as fit_log:
            log_lines = fit_log.readlines()

        rmsd = float(log_lines[-6].split()[2])

        system.format_print("Completed fit number {} with rmsd {}.".format(attempts, rmsd), italics=True)

        system.format_print("Current best fit has rmsd {}.".format(best_rmsd), italics=True)

        # if the new fit is better than the old fit, replace the best log and best cdl files
        if rmsd < best_rmsd:
            system.format_print("Replaced previous best fit with most recent one.", italics=True)

            os.rename(fit_log_path, best_fit_log_path)
            os.rename("fit-1b.cdl", "best-fit-1b.cdl")
            os.rename("fit-1b-initial.cdl", "best-fit-1b-initial.cdl")
            os.rename("correlation.dat", "best-correlation.dat")

            best_rmsd = rmsd
            
        attempts += 1

        system.format_print("\n", italics=True)

    # remove the most recent fit file
    try:
        os.remove(fit_log_path)
        os.remove("fit-1b.cdl")
        os.remove("fit-1b-initial.cdl")
        os.remove("correlation.dat")
    # in the case that there is no most recent fit file because the last fit was the best fit, do nothing
    except FileNotFoundError:
        pass

    system.format_print("Completed {} fits.".format(num_fits), bold=True, color=system.Color.GREEN)


def create_1b_nc_file(settings_path, fit_dir_path, fitted_nc_path):
    
    system.call("ncgen", "-o", fitted_nc_path, "best-fit-1b.cdl")

    os.rename("best-fit-1b.cdl", os.path.join(fit_dir_path, "fit-1b.cdl"))
    os.rename("best-fit-1b-initial.cdl", os.path.join(fit_dir_path, "fit-1b-initial.cdl"))
    os.rename("best-correlation.dat", os.path.join(fit_dir_path, "correlation.dat"))


def fit_1b_training_set(settings_path, fit_code_path, training_set_path, fit_dir_path, fitted_nc_path, num_fits = 10):
    """
    Fits a given 1b training set using a given 1b fit code.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fit_code_path       - Local path to the fit code with which to fit the training set. Generated in fit_dir_path
                by compile_fit_code.
        training_set_path   - Local path to the file to read the training set from.
        fit_dir_path        - Local path to the directory to write the ".cdl" and ".dat" files created by this fit.
        fitted_nc_path      - Local path to the file to write the final fitted ".nc" to.
        num_fits            - Performs this many fits and then gives only the best one.

    Returns:
        None
    """

    perform_1b_fits(settings_path, fit_code_path, training_set_path, fit_dir_path, num_fits = num_fits)
    create_1b_nc_file(settings_path, fit_dir_path, fitted_nc_path)


def perform_2b_ttm_fits(settings_path, fit_code_path, training_set_path, fit_dir_path, num_fits = 10):

    system.format_print("Performing {} fits from which the best will be chosen...".format(num_fits), bold=True, color=system.Color.GREEN)

    settings = SettingsReader(settings_path)

    files.init_directory(fit_dir_path)

    best_fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "ttm", "best_fit.log"))
    fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "ttm", "fit.log"))

    attempts = 1
    with open(best_fit_log_path, "w") as best_fit_log:
        system.call(fit_code_path, training_set_path, out_file = best_fit_log)
        os.rename("individual_terms.dat", "best-individual_terms.dat")
        os.rename("ttm-params.txt", "best-ttm-params.txt")
        os.rename("correlation.dat", "best-correlation.dat")

    with open(best_fit_log_path, "r") as best_fit_log:
        best_log_lines = best_fit_log.readlines()

    best_rmsd = float(best_log_lines[-4].split()[2])

    system.format_print("Completed first fit with rmsd {}.\n".format(best_rmsd), italics=True)

    while(attempts < num_fits):

        with open(fit_log_path, "w") as fit_log:
            system.call(fit_code_path, training_set_path, out_file = fit_log)
        
        with open(fit_log_path, "r") as fit_log:
            log_lines = fit_log.readlines()

        rmsd = float(log_lines[-4].split()[2])

        system.format_print("Completed fit number {} with rmsd {}.".format(attempts, rmsd), italics=True)

        system.format_print("Current best fit has rmsd {}.".format(best_rmsd), italics=True)

        if rmsd < best_rmsd:

            system.format_print("Replaced previous best fit with most recent one.")

            os.rename(fit_log_path, best_fit_log_path)
            os.rename("individual_terms.dat", "best-individual_terms.dat")
            os.rename("ttm-params.txt", "best-ttm-params.txt")
            os.rename("correlation.dat", "best-correlation.dat")

            best_rmsd = rmsd
            

        attempts += 1

        system.format_print("\n", italics=True)

    # remove the most recent fit file
    try:
        os.remove(fit_log_path)
        os.remove("individual_terms.dat")
        os.remove("ttm-params.txt")
        os.remove("correlation.dat")
    # in the case that there is no most recent fit file because the last fit was the best fit, do nothing
    except FileNotFoundError:
        pass

    os.rename("best-individual_terms.dat", os.path.join(fit_dir_path, "individual_terms.dat"))
    os.rename("best-ttm-params.txt", os.path.join(fit_dir_path, "ttm-params.txt"))
    os.rename("best-correlation.dat", os.path.join(fit_dir_path, "correlation.dat"))

    system.format_print("Completed {} fits.".format(num_fits), bold=True, color=system.Color.GREEN)

def add_A_and_b_to_config_file(settings_path, fit_dir_path, config_path):
    
    # read d6 and A constants from ttm output file
    with open(os.path.join(fit_dir_path, "ttm-params.txt"), "r") as ttm_file:
        A = [float(a) for a in ttm_file.readline().split()]
        d6 = [float(d) for d in ttm_file.readline().split()]

    # write d6 and A to the config.ini file
    lines = []
    with open(config_path, "r") as config_file:
        for line in config_file:
            lines.append(line)

    found_A = False
    with open(config_path, "w") as config_file:
        for line in lines:
            if line.startswith("A = "):
                found_A = True
            if not line.startswith("d6 = "):
                config_file.write(line)

            else:
                config_file.write("d6 = {}\n".format([[], [], d6]))

        if not found_A:
            config_file.write("A = {}\n".format([[], [], A]))

def fit_2b_ttm_training_set(settings_path, fit_code_path, training_set_path, fit_dir_path, config_path, num_fits = 10):
    """
    Fits a given 2b training set using a given 2b ttm fit code

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fit_code_path       - Local path to the fit code with which to fit the training set. Generated in fit_dir_path
                by compile_fit_code.
        training_set_path   - Local path to the file to read the training set from.
        fit_dir_path        - the directory where the fit log and other files created by the fit go
        config_path         - Local path to the ".ini" file to write A6 and d6 constants to.
        num_fits            - Performs this many fits and then gives only the best one.

    Returns:
        None
    """

    perform_2b_ttm_fits(settings_path, fit_code_path, training_set_path, fit_dir_path, num_fits = num_fits)
    add_A_and_b_to_config_file(settings_path, fit_dir_path, config_path)


def perform_2b_fits(settings_path, fit_code_path, training_set_path, fit_dir_path, num_fits = 10):

    system.format_print("Performing {} fits from which the best will be chosen...".format(num_fits), bold=True, color=system.Color.YELLOW)
    
    settings = SettingsReader(settings_path)

    # init the required files
    files.init_directory(fit_dir_path)
    best_fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "mb", "best-fit.log"))
    fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "mb" "fit.log"))

    # keeps tracks of the number of attempts to generate a fit
    attempts = 1

    # generate an initial fit
    with open(best_fit_log_path, "w") as best_fit_log:
        system.call(fit_code_path, training_set_path, out_file = best_fit_log)
        os.rename("fit-2b.cdl", "best-fit-2b.cdl")
        os.rename("fit-2b-initial.cdl", "best-fit-2b-initial.cdl")
        os.rename("correlation.dat", "best-correlation.dat")

    with open(best_fit_log_path, "r") as best_fit_log:
        best_log_lines = best_fit_log.readlines()

    best_rmsd = float(best_log_lines[-6].split()[2])

    system.format_print("Completed first fit with rmsd {}.\n".format(best_rmsd), italics=True)

    while(attempts < num_fits):

        # generate a new fit
        with open(fit_log_path, "w") as fit_log:
            system.call(fit_code_path, training_set_path, out_file = fit_log)
  
        with open(fit_log_path, "r") as fit_log:
            log_lines = fit_log.readlines()

        rmsd = float(log_lines[-6].split()[2])

        system.format_print("Completed fit number {} with rmsd {}.".format(attempts, rmsd), italics=True)

        system.format_print("Current best fit has rmsd {}.".format(best_rmsd), italics=True)

        # if the new fit is better than the old fit, replace the best log and best cdl files
        if rmsd < best_rmsd:
            system.format_print("Replaced previous best fit with most recent one.", italics=True)

            os.rename(fit_log_path, best_fit_log_path)
            os.rename("fit-2b.cdl", "best-fit-2b.cdl")
            os.rename("fit-2b-initial.cdl", "best-fit-2b-initial.cdl")
            os.rename("correlation.dat", "best-correlation.dat")

            best_rmsd = rmsd
            
        attempts += 1

        system.format_print("\n", italics=True)

    # remove the most recent fit file
    try:
        os.remove(fit_log_path)
        os.remove("fit-2b.cdl")
        os.remove("fit-2b-initial.cdl")
        os.remove("correlation.dat")
    # in the case that there is no most recent fit file because the last fit was the best fit, do nothing
    except FileNotFoundError:
        pass

    system.format_print("Completed {} fits.".format(num_fits), bold=True, color=system.Color.GREEN)

def create_2b_nc_file(settings_path, fit_dir_path, fitted_nc_path):
    
    system.call("ncgen", "-o", fitted_nc_path, "best-fit-2b.cdl")

    os.rename("best-fit-2b.cdl", os.path.join(fit_dir_path, "fit-2b.cdl"))
    os.rename("best-fit-2b-initial.cdl", os.path.join(fit_dir_path, "fit-2b-initial.cdl"))
    os.rename("best-correlation.dat", os.path.join(fit_dir_path, "correlation.dat"))

def fit_2b_training_set(settings_path, fit_code_path, training_set_path, fit_dir_path, fitted_nc_path, num_fits = 10):
    """
    Fits the fit code to a given training set

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fit_code            - Local path to the code to fit.
        training_set        - Local path to the training set to fit the code to.
        fit_dir_path        - Local path to the directory where the .cdl and .dat files will end up.
        fitted_nc_path      - Local path to file to write final fitted ".nc" file to.
        num_fits            - Performs this many fits and then gives only the best one.

    Returns:
        None
    """

    perform_2b_fits(settings_path, fit_code_path, training_set_path, fit_dir_path, num_fits = num_fits)
    create_2b_nc_file(settings_path, fit_dir_path, fitted_nc_path)

def prepare_fits(settings_path, fit_executable_path, training_set_path, DE = 20, alpha = 0.0005, num_fits = 10, ttm = False):
    # Get information
    settings = SettingsReader(settings_path)
    workdir = os.getcwd()
    if ttm:
        fit_folder_name = "ttm_nrg_fits"
    else:
        fit_folder_name = "mb_nrg_fits"
    fit_folder_prefix = workdir + "/" + settings.get("files", "log_path") + "/" + fit_folder_name + "/"
    if not os.path.exists(fit_folder_prefix):
        os.mkdir(fit_folder_prefix)
    
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
            # Link the training set to the fit folder
            system.call("ln", "-s", "{}/{}".format(workdir, training_set_path), ".")
            # Create bash script that will run the fit
            my_bash = open("run_fit.sh",'w')
            my_bash.write("#!/bin/bash\n")
            my_bash.write("\n{}/{} {} {} {} > fit.log 2> fit.err \n".format(workdir,fit_executable_path,training_set_path, DE, alpha))  
            my_bash.close()
            system.call("chmod", "744", "run_fit.sh")
            fit_index += 1
            count += 1 
            print("Succesfully created fit folder {}.".format(new_fit_folder))

    os.chdir(workdir)

def execute_fits(settings_path, ttm = False):
    # Get information
    settings = SettingsReader(settings_path)
    workdir = os.getcwd()
    if ttm:
        fit_folder_name = "ttm_nrg_fits"
    else:
        fit_folder_name = "mb_nrg_fits"
    fit_folder_prefix = workdir + "/" + settings.get("files", "log_path") + "/" + fit_folder_name + "/"

    # A fit folder that does not have a fit.log file inside is considered not run
    # In that case, the run_fit.sh will be executed
    os.chdir(fit_folder_prefix)
    all_fits = glob.glob("fit*")
    for fit in all_fits:
        os.chdir(fit)
        if not os.path.exists("fit.log"):
            print("{} is running.".format(fit))
            system.call("./run_fit.sh")
            print("{} is completed.".format(fit))
        else:
            print("{} is already done. Continuing...".format(fit))
        os.chdir("../")

    os.chdir(workdir)
    
def retrieve_best_fit(settings_path, ttm = False, fitted_nc_path = "mbnrg.nc"):
    # Get information
    settings = SettingsReader(settings_path)
    workdir = os.getcwd()
    if ttm:
        fit_folder_name = "ttm_nrg_fits"
    else:
        fit_folder_name = "mb_nrg_fits"
    fit_folder_prefix = workdir + "/" + settings.get("files", "log_path") + "/" + fit_folder_name + "/"

    # Loop over all the fits, check the output, and store the results.
    os.chdir(fit_folder_prefix)
    all_fits = glob.glob("fit*")
    results = []
    for fit in all_fits:
        os.chdir(fit)

        try:
            with open("fit.log",'r') as logfile:
                log_lines = logfile.readlines()
                full_rmsd = float(log_lines[-7].split()[2])
                wfull_rmsd = float(log_lines[-6].split()[2])
                max_error = float(log_lines[-5].split()[2])
                low_rmsd = float(log_lines[-4].split()[2])
                low_max_error = float(log_lines[-3].split()[2])
                results.append([fit, full_rmsd, wfull_rmsd, max_error, low_rmsd, low_max_error])
        except:
            print("Doesn't seem that the log file in " + fit + " is correct...")
            print("Maybe you want to rerun " + fit + " again.")
            results.append([fit, float('inf'), float('inf'), float('inf'), float('inf'), float('inf')])

        os.chdir("../")

    # Sort the results according to weighet RMSD of the full TS
    sorted_results = sorted(results, key=lambda x: x[2])

    # Get the best result
    best_fit = sorted_results[0][0]
    best_results = sorted_results[0]
    
    # Check if best fit folder exists.
    # If it does exist, check the output in best folder to ensure that the new best
    # fit is actually the best
    # If not, just store the best fit in there
    if not os.path.exists("best_fit"):
        system.call("cp", "-r", best_fit, "best_fit")
    else:
        os.chdir("best_fit")
        with open("fit.log",'r') as logfile:
            log_lines = logfile.readlines()
            full_rmsd = float(log_lines[-7].split()[2])
            wfull_rmsd = float(log_lines[-6].split()[2])
            max_error = float(log_lines[-5].split()[2])
            low_rmsd = float(log_lines[-4].split()[2])
            low_max_error = float(log_lines[-3].split()[2])
            previous_best = [fit, full_rmsd, wfull_rmsd, max_error, low_rmsd, low_max_error]
        os.chdir("../")

        if previous_best[2] <= sorted_results[0][2]:
            print("Previous fit is better than any fit tested")
            best_results = previous_best
        else:
            print("Replacing best_fit by {}".format(sorted_results[0][0]))
            system.call("rm", "-rf", "best_fit")
            system.call("cp", "-r", sorted_results[0][0], "best_fit")

    os.chdir("best_fit")
    nb = len(settings.get("molecule","names").split(","))
    if os.path.exists("fit-" + str(nb) + "b.cdl"):
        system.call("ncgen", "-o", fitted_nc_path, "fit-" + str(nb) + "b.cdl")

    # Report best RMSD
    print("Best fit found has a weighted RMSD of {} kcal/mol, a low energy RMSD of {} kcal/mol, and a maximum error in the low energy training set of {} kcal/mol".format(best_results[2], best_results[4], best_results[5]))

    os.chdir(workdir)

def update_config_with_ttm(settings_path, config_path):
    config = SettingsReader(config_path) 
    settings = SettingsReader(settings_path)
    workdir = os.getcwd()
    fit_folder_name = "ttm_nrg_fits"
    fit_folder_prefix = workdir + "/" + settings.get("files", "log_path") + "/" + fit_folder_name + "/best_fit/"

    with open(fit_folder_prefix + "ttm-nrg_params.dat",'r') as ttm_file:
        a_buck = ttm_file.readline().strip().split()
        b_buck = ttm_file.readline().strip().split()

    config.set("fitting","A",a_buck)    
    config.set("fitting","d6",b_buck)   

    config.write(config_path)












