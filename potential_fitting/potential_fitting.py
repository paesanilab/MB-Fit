import os, subprocess, contextlib
import random

from potential_fitting.utils import SettingsReader
from . import configurations, database, polynomials, fitting
from .database import Database

def create_dirs(settings_path):
    """
    Creates the log directory

    Args:
        settings_path   - the file containing all relevent settings information

    Returns:
        None
    """
    settings = SettingsReader(settings_path)

    if not os.path.isdir(settings.get("files", "log_path")):
        os.mkdir(settings.get("files", "log_path"))

def optimize_geometry(settings_path, unopt_geo, opt_geo):
    """
    Optimizes the geometry of the given molecule

    Args:
        settings_path    - the file containing all relevent settings information
        unopt_geo   - file to read the unoptimized geoemtry
        opt_geo     - file to write the optimized geometry

    Returns:
        None
    """

    if os.path.dirname(opt_geo) == "":
        opt_geo = "./" + opt_geo
    if not os.path.isdir(os.path.dirname(opt_geo)):
        os.mkdir(os.path.dirname(opt_geo))

    configurations.optimize_geometry(settings_path, unopt_geo, opt_geo)

def generate_normal_modes(settings_path, geo, normal_modes):
    """
    Generates the normal modes for the given molecule

    Args:
        settings_path - the file containing all relevent settings information
        geo         - file to read optimized geometry
        normal_modes - file to write normal modes

    Returns:
        dim_null of normal modes

    Other Effects:
        Writes DIM NULL to console
    """
    
    if os.path.dirname(normal_modes) == "":
        normal_modes = "./" + normal_modes
    if not os.path.isdir(os.path.dirname(normal_modes)):
        os.mkdir(os.path.dirname(normal_modes))

    dim_null = configurations.generate_normal_modes(settings_path, geo, normal_modes)

    return dim_null


def generate_1b_configurations(settings_path, geo, normal_modes, dim_null, config_path):
    """
    Generates 1b configurations for a given monomer from a set of normal modes

    Args:
        settings_path - the file containing all relevent settings information
        geo         - file to read optimized geometry
        normal_modes - file to read normal modes
        dim_null    - the DIM NULL of this molecule, see generate_normal_modes()
        config_path - file to write configurations

    Returns:
        None
    """

    if os.path.dirname(config_path) == "":
        config_path = "./" + config_path
    if not os.path.isdir(os.path.dirname(config_path)):
        os.mkdir(os.path.dirname(config_path))

    configurations.generate_1b_configurations(settings_path, geo, normal_modes, dim_null, config_path)

def generate_2b_configurations(settings_path, geo1, geo2, number_of_configs, config_path, min_distance = 1, max_distance = 5, min_inter_distance = 1.2, use_grid = False, step_size = 0.5, seed = random.randint(-1000000, 1000000)):
    """
    Generates 2b configurations for a given dimer

    Args:
        settings_path - the file containing all relevent settings information
        geo1        - the first optimized geometry
        geo2        - the second optimized geometry
        number_of_configs - target number of configurations to generate, if max_distance is set too low or min_inter_distance is set too high, then less configurations may be generated
        config_path - path to file in which to generate configurations
        min_distance - minimum distance between the centers of mass of the two molecules
        max_distance - the maximum distance between the centers of mass of the two molecules
        min_inter_distance - the minimum distance of any intermolecular pair of atoms
        use_grid - if False, configurations are space roughly evenly between min_distance and max_distance. If True, then configurations are placed at intervals along this distance.
        step_size - if use_grid is True, then this dictates the distance of the spacing interval used to place the centers of masses of the molecules, otherwise, this parameter has no effect.
        seed - the same seed will generate the same configurations.

    Returns:
        None
    """

    if os.path.dirname(config_path) == "":
        config_path = "./" + config_path
    if not os.path.isdir(os.path.dirname(config_path)):
        os.mkdir(os.path.dirname(config_path))

    configurations.generate_2b_configurations(geo1, geo2, number_of_configs, config_path, min_distance, max_distance, min_inter_distance, use_grid, step_size, seed)

def init_database(settings_path, database_name, config_files):
    """
    Creates a database from the given config files. Can be called on a new file
    to create a new database, or an existing database to add more energies to be calculated

    Args:
        settings_path - the file containing all relevent settings information
        database_name - file to make this database in
        config_files - directory with .xyz and optimized geometry .opt.xyz inside

    Returns:
        None
    """

    if os.path.dirname(database_name) == "":
        database_name = "./" + database_name
    if not os.path.isdir(os.path.dirname(database_name)):
        os.mkdir(os.path.dirname(database_name))

    database.initialize_database(settings_path, database_name, config_files)

def fill_database(settings_path, database_name):
    """
    Fills a given database with calculated energies, MAY TAKE A WHILE
    
    Args:
        settings_path - the file containing all relevent settings information
        database_name - the file in which the database is stored

    Returns:
        None
    """

    database.fill_database(settings_path, database_name)

def generate_1b_training_set(settings_path, database_name, training_set, molecule_name, method = "%", basis = "%", cp = "%", tag = "%"):
    """
    Generates a training set from the energies inside a database.

    Specific method, bnasis, and cp may be specified to only use energies calculated
    with a specific model.

    '%' can be used to stand in as a wild card, meaning any method/basis/cp is ok

    Args:
        settings_path    - the file containing all relevent settings information
        database_name - the file in which the database is stored
        training_set - the file to write the training set to
        molecule_name - the name of the moelcule to generate a training set for
        method      - only use energies calcualted by this method
        basis       - only use energies calculated in this basis
        cp          - only use energies calculated with the same cp
        tag         - only use energies marked with this tag

    Returns:
        None
    """

    if os.path.dirname(training_set) == "":
        training_set = "./" + training_set
    if not os.path.isdir(os.path.dirname(training_set)):
        os.mkdir(os.path.dirname(training_set))

    database.generate_1b_training_set(settings_path, database_name, training_set, molecule_name, method, basis, cp, tag)

def generate_2b_training_set(settings_path, database_name, training_set, monomer1_name, monomer2_name, method = "%", basis = "%", cp = "%", tag = "%"):
    """
    Generates a training set from the energies inside a database.

    Specific method, bnasis, and cp may be specified to only use energies calculated
    with a specific model.

    '%' can be used to stand in as a wild card, meaning any method/basis/cp is ok

    Args:
        settings_path    - the file containing all relevent settings information
        database_name - the file in which the database is stored
        training_set - the file to write the training set to
        monomer1_name - name of first monomer in the dimer
        monomer2_name - name of the second monomer in the dimer
        method      - only use energies calcualted by this method
        basis       - only use energies calculated in this basis
        cp          - only use energies calculated with the same cp
        tag         - only use energies marked with this tag

    Returns:
        None
    """

    if os.path.dirname(training_set) == "":
        training_set = "./" + training_set
    if not os.path.isdir(os.path.dirname(training_set)):
        os.mkdir(os.path.dirname(training_set))

    database.generate_2b_training_set(settings_path, database_name, training_set, monomer1_name, monomer2_name, method, basis, cp, tag)

def generate_poly_input(settings_path, poly_in_path):
    """
    Generates an input file for the polynomial generation script

    Args:
        settings_path    - the file containing all relevent settings information
        poly_in_path - the file to write the polynomial generation input
                     name of file should be in format A1B2.in, it is ok to have
                     extra directories prior to file name (thisplace/thatplace/A3.in)

    Returns:
        None
    """

    if os.path.dirname(poly_in_path) == "":
        poly_in_path = "./" + poly_in_path
    if not os.path.isdir(os.path.dirname(poly_in_path)):
        os.mkdir(os.path.dirname(poly_in_path))

    polynomials.generate_input_poly(settings_path, poly_in_path)

def generate_poly_input_from_database(settings_path, database_name, molecule_name, poly_directory):
    """
    Generates an input file for the polynomial generation script.
    Looks in a database to find the symmetry and creates a file in the given directory.

    Args:
        settings_path - the file containing all relevent settings information
        database_name - the name of the database to read the symmetry from
        molecule_name - name of molecule to read symmetry of
        poly_directory - the directory to make the input file in

    Returns:
        None
    """

    with Database(database_name) as database:
        symmetry = database.get_symmetry(molecule_name)

        generate_poly_input(settings_path, poly_directory + "/" + symmetry + ".in")

def generate_polynomials(settings_path, poly_in_path, order, poly_directory):
    """
    Generates the maple and cpp polynomial codes

    Args:
        settings_path    - the file containing all relevent settings information
        poly_in_path - the file to read polynomial input from
        order - the order of the polynomial to generate
        poly_directory - the directory to place the polynomial files in

    Returns:
        None
    """

    this_file_path = os.path.dirname(os.path.abspath(__file__))

    original_dir = os.getcwd()

    if not os.path.isdir(poly_directory):
        os.mkdir(poly_directory)

    os.chdir(poly_directory)

    os.system("poly-gen_mb-nrg.pl " + str(order) + " " + original_dir + "/" + poly_in_path + " > poly.log")

    os.chdir(original_dir)

def execute_maple(settings_path, poly_directory):
    """
    Runs maple on the polynomial files to turn them into actual cpp files

    Args:
        settings_path    - the file containing all relevent settings information
        poly_directory - the directory with all the polynomial files
    """

    this_file_path = os.path.dirname(os.path.abspath(__file__))

    original_dir = os.getcwd()

    os.chdir(poly_directory)

    # clear any old files, because for some reason maple appends to existing files instead of clobbering
    with contextlib.suppress(FileNotFoundError):
        os.remove("poly-grd.c")
        os.remove("poly-nogrd.c")
        os.remove("poly-grd.cpp")
        os.remove("poly-nogrd.cpp")
    
    os.system("maple poly-grd.maple")
    os.system("maple poly-nogrd.maple")

    os.system(this_file_path + "/../polynomial_generation/src/clean-maple-c.pl < poly-grd.c > poly-grd.cpp")
    os.system(this_file_path + "/../polynomial_generation/src/clean-maple-c.pl < poly-nogrd.c > poly-nogrd.cpp")

    os.chdir(original_dir)

def generate_fit_config(settings_path, molecule_in, config_path, *opt_geometry_paths, distance_between = 20):
    """
    Generates the fit config for the optimized geometry

    For >=2b, the optimized geometry should be the optimized geometry of each monomer seperated by roughly
    40 angstroms

    Args:
        settings_path   - the file containing all relevent settings information
        molecule_in     - string of fromat "A1B2"
        config_path     - path to file to write fit config in
        opt_geometry_paths - the paths to the optimized geometries to include in this fit config, should be 1 to 3 (inclusive)
        distance_between - distance between each geometry in the qchem calculation. If the qchemc calculation does not converge, try different values of this

    Returns:
        None
    """

    if os.path.dirname(config_path) == "":
        config_path = "./" + config_path
    if not os.path.isdir(os.path.dirname(config_path)):
        os.mkdir(os.path.dirname(config_path))

    fitting.make_config(settings_path, molecule_in, config_path, *opt_geometry_paths, distance_between = distance_between)
    
def generate_1b_fit_code(settings_path, config, poly_in_path, poly_path, poly_order, fit_directory):
    """
    Generates the fit code based on the polynomials for a monomer

    Args:
        settings_path - the file containing all relevent settings information
        config    - monomer config file
        poly_in_path - the A3B2.in type file
        poly_path   - directory where polynomial files are
        poly_order - the order of the polynomial in poly_path
        fit_directory - directory to generate fit code in

    Returns:
        None
    """

    if not os.path.isdir(fit_directory):
        os.mkdir(fit_directory)
    
    fitting.prepare_1b_fitting_code(config, poly_in_path, poly_path, poly_order, fit_directory)

def generate_2b_ttm_fit_code(settings_path, config, molecule_in, fit_directory):
    """
    Generates the fit TTM fit code for a dimer

    Args:
        settings_path - the file containing all relevent settings information
        config    - config file generated by generate_fit_config (requires manual editing)
        poly_in_path - the A3B2.in name for this dimer
        fit_directory - directory to generate fit code in

    Returns:
        None
    """

    this_file_path = os.path.dirname(os.path.abspath(__file__))

    original_dir = os.getcwd()

    if not os.path.isdir(fit_directory):
        os.mkdir(fit_directory)

    os.chdir(fit_directory) 
    
    subprocess.call("cp " + this_file_path + "/../codes/2b-codes/template/* .", shell=True)
    subprocess.call("python " + this_file_path + "/../codes/2b-codes/get_2b_TTM_codes.py {} {} {}".format(original_dir + "/" + settings_path, original_dir + "/" + config, molecule_in), shell=True)

    os.chdir(original_dir)   


def compile_fit_code(settings_path, fit_directory):
    """
    Compiles the fit code in the given directory

    Args:
        settings_path    - the file containing all relevent settings information
        fit_directory - the directory with the fit code

    Returns:
        None
    """

    original_dir = os.getcwd()

    os.chdir(fit_directory)
    os.system("make clean")
    os.system("make")
    os.chdir(original_dir)

def fit_1b_training_set(settings_path, fit_code, training_set, fit_directory, fitted_code):
    """
    Fits the fit code to a given training set

    Args:
        settings_path    - the file containing all relevent settings information
        fit_code    - the code to fit
        training_set - the training set to fit the code to
        fit_directory - the directory where the .cdl and .dat files will end up
        fitted_code - file to write final fitted code to

    Returns:
        None
    """

    if not os.path.isdir(fit_directory):
        os.mkdir(fit_directory)
    if os.path.dirname(fitted_code) == "":
        fitted_code = "./" + fitted_code
    if not os.path.isdir(os.path.dirname(fitted_code)):
        os.mkdir(os.path.dirname(fitted_code))

    os.system(fit_code + " " + training_set)

    os.system("ncgen -o fit-1b.nc fit-1b.cdl")

    os.rename("fit-1b.cdl", os.path.join(fit_directory, "fit-1b.cdl"))
    os.rename("fit-1b-initial.cdl", os.path.join(fit_directory, "fit-1b-initial.cdl"))
    os.rename("correlation.dat", os.path.join(fit_directory, "correlation.dat"))
    os.rename("fit-1b.nc", fitted_code)

def fit_2b_ttm_training_set(settings_path, fit_code, training_set, fit_directory):
    """
    Fits the ttm fit code to a given training set

    Args:
        settings_path    - the file containing all relevent settings information
        fit_code    - the code to fit
        training_set - the training set to fit the code to
        fit_directory - the directory where the fit log and other files created by the fit go

    Returns:
        None
    """

    if not os.path.isdir(fit_directory):
        os.mkdir(fit_directory)

    attempts = 1;
    os.system(fit_code + " " + training_set + " > " + fit_directory + "/best_fit.log")
    while(attempts < 10):
        os.system(fit_code + " " + training_set + " > " + fit_directory + "/fit.log")
        
        with open(fit_directory + "/fit.log", "r") as fit_log, open(fit_directory + "/best_fit.log", "r") as best_fit_log:
            log_lines = fit_log.readlines()
            best_log_lines = best_fit_log.readlines()

        rmsd = float(log_lines[-4].split()[2])
        best_rmsd = float(best_log_lines[-4].split()[2])

        if rmsd < best_rmsd:
            os.rename(os.path.join(fit_directory, "fit.log"), os.path.join(fit_directory, "best_fit.log"))
            

        attempts += 1

    os.remove(os.path.join(fit_directory, "fit.log"))
    os.rename("individual_terms.dat", os.path.join(fit_directory, "individual_terms.dat"))
    os.rename("ttm-params.txt", os.path.join(fit_directory, "ttm-params.txt"))
    os.rename("correlation.dat", os.path.join(fit_directory, "correlation.dat"))
