import sys, os, subprocess
import random, configparser

def create_dirs(settings):
    """
    Creates the log directory

    Args:
        settings    - the file containing all relevent settings information

    Returns:
        None
    """

    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)

    if not os.path.isdir(config['files']['log_path']):
        os.mkdir(config['files']['log_path'])

def optimize_geometry(settings, unopt_geo, opt_geo):
    """
    Optimizes the geometry of the given molecule

    Args:
        settings    - the file containing all relevent settings information
        unopt_geo   - file to read the unoptimized geoemtry
        opt_geo     - file to write the optimized geometry

    Returns:
        None
    """

    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../configuration_generator/src")
    import geometry_optimizer
   
    if os.path.dirname(opt_geo) == "":
        opt_geo = "./" + opt_geo
    if not os.path.isdir(os.path.dirname(opt_geo)):
        os.mkdir(os.path.dirname(opt_geo))

    geometry_optimizer.optimize_geometry(settings, unopt_geo, opt_geo)

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../configuration_generator/src") 

def generate_normal_modes(settings, geo, normal_modes):
    """
    Generates the normal modes for the given molecule

    Args:
        settigns    - the file containing all relevent settings information
        geo         - file to read optimized geometry
        normal_modes - file to write normal modes

    Returns:
        dim_null of normal modes

    Other Effects:
        Writes DIM NULL to console
    """
    
    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../configuration_generator/src")
    import normal_modes_generator
   
    if os.path.dirname(normal_modes) == "":
        normal_modes = "./" + normal_modes
    if not os.path.isdir(os.path.dirname(normal_modes)):
        os.mkdir(os.path.dirname(normal_modes))

    dim_null = normal_modes_generator.generate_normal_modes(settings, geo, normal_modes)

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../configuration_generator/src") 

    return dim_null


def generate_1b_configurations(settings, geo, normal_modes, dim_null, configurations):
    """
    Generates 1b configurations for a given monomer from a set of normal modes

    Args:
        settings    - the file containing all relevent settings information
        geo         - file to read optimized geometry
        normal_modes - file to read normal modes
        dim_null    - the DIM NULL of this molecule, see generate_normal_modes()
        configurations - file to write configurations

    Returns:
        None
    """

    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../configuration_generator/src")
    import configuration_generator
   
    if not os.path.isdir(configurations):
        os.mkdir(configurations)

    configuration_generator.generate_configurations(settings, geo, normal_modes, dim_null, configurations)

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../configuration_generator/src") 

def generate_2b_configurations(settings, geo1, geo2, configs_per_distance, min_distance, configurations, flex = False, seed = random.randint(-1000000, 1000000)):
    """
    Generates 2b configurations for a given dimer

    Args:
        settings    - the file containing all relevent settings information
        geo1        - the first optimized geometry
        geo2        - the second optimized geometry
        configs_per_distance - number of configurations per distance
        min_distance - min distance between any 2 molecules from each monomer
        configurations - file to write configurations
        flex        - True will generate flex rather than rigid configurations. Default is False
        seed        - seed to generate random configs, the same seed will yeild the same configurations. Defualt is a random seed

    Returns:
        None
    """

    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)

    original_dir = os.getcwd()

    if not os.path.isdir(configurations):
        os.mkdir(configurations)

    os.chdir(configurations)

    if flex:
        os.system(os.path.dirname(os.path.abspath(__file__)) + "/../configuration_generator/2b_configs/bin/2b_ts_flex " + original_dir + "/" + geo1 + " " + original_dir + "/" + geo2 + " " + str(configs_per_distance) + " " + str(min_distance) + " " + str(seed) + " > /dev/null")
    else:
        os.system(os.path.dirname(os.path.abspath(__file__)) + "/../configuration_generator/2b_configs/bin/2b_ts_rigid " + original_dir + "/" + geo1 + " " + original_dir + "/" + geo2 + " " + str(configs_per_distance) + " " + str(min_distance) + " " + str(seed) + " > /dev/null")

    os.chdir(original_dir)


def init_database(settings, database_name, config_files):
    """
    Creates a database from the given config files. Can be called on a new file
    to create a new database, or an existing database to add more energies to be calculated

    Args:
        settings    - the file containing all relevent settings information
        database_name - file to make this database in
        config_files - directory with .xyz and optimized geometry .opt.xyz inside

    Returns:
        None
    """

    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../qm_mb_energy_calculator/src")
    import database_initializer

    if os.path.dirname(database_name) == "":
        database_name = "./" + database_name
    if not os.path.isdir(os.path.dirname(database_name)):
        os.mkdir(os.path.dirname(database_name))

    database_initializer.initialize_database(settings, database_name, config_files)

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../qm_mb_energy_calculator/src") 

def fill_database(settings, database_name):
    """
    Fills a given database with calculated energies, MAY TAKE A WHILE
    
    Args:
        settings    - the file containing all relevent settings information
        database_name - the file in which the database is stored

    Returns:
        None
    """

    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../qm_mb_energy_calculator/src")
    import database_filler

    database_filler.fill_database(settings, database_name)

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../qm_mb_energy_calculator/src") 

def generate_1b_training_set(settings, database_name, training_set, molecule_name, method = "%", basis = "%", cp = "%", tag = "%"):
    """
    Generates a training set from the energies inside a database.

    Specific method, bnasis, and cp may be specified to only use energies calculated
    with a specific model.

    '%' can be used to stand in as a wild card, meaning any method/basis/cp is ok

    Args:
        settings    - the file containing all relevent settings information
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

    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../qm_mb_energy_calculator/src")
    import training_set_generator

    if os.path.dirname(training_set) == "":
        training_set = "./" + training_set
    if not os.path.isdir(os.path.dirname(training_set)):
        os.mkdir(os.path.dirname(training_set))

    training_set_generator.generate_1b_training_set(settings, database_name, training_set, molecule_name, method, basis, cp, tag)

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../qm_mb_energy_calculator/src") 
 
def generate_2b_training_set(settings, database_name, training_set, monomer1_name, monomer2_name, method = "%", basis = "%", cp = "%", tag = "%"):
    """
    Generates a training set from the energies inside a database.

    Specific method, bnasis, and cp may be specified to only use energies calculated
    with a specific model.

    '%' can be used to stand in as a wild card, meaning any method/basis/cp is ok

    Args:
        settings    - the file containing all relevent settings information
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

    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../qm_mb_energy_calculator/src")
    import training_set_generator

    if os.path.dirname(training_set) == "":
        training_set = "./" + training_set
    if not os.path.isdir(os.path.dirname(training_set)):
        os.mkdir(os.path.dirname(training_set))

    training_set_generator.generate_2b_training_set(settings, database_name, training_set, monomer1_name, monomer2_name, method, basis, cp, tag)

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../qm_mb_energy_calculator/src") 
    


def generate_poly_input(settings, poly_in_path):
    """
    Generates an input file for the polynomial generation script

    Args:
        settings    - the file containing all relevent settings information
        poly_in_path - the file to write the polynomial generation input
                     name of file should be in format A1B2.in, it is ok to have
                     extra directories prior to file name (thisplace/thatplace/A3.in)

    Returns:
        None
    """

    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../polynomial_generation/src")
    import generate_input_poly


    if os.path.dirname(poly_in_path) == "":
        poly_in_path = "./" + poly_in_path
    if not os.path.isdir(os.path.dirname(poly_in_path)):
        os.mkdir(os.path.dirname(poly_in_path))

    generate_input_poly.generate_input_poly(settings, poly_in_path)

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../polynomial_generation/src") 


def generate_polynomials(settings, poly_in_path, order, poly_directory):
    """
    Generates the maple and cpp polynomial codes

    Args:
        settings    - the file containing all relevent settings information
        poly_in_path - the file to read polynomial input from
        order - the order of the polynomial to generate
        poly_directory - the directory to place the polynomial files in

    Returns:
        None
    """

    this_file_path = os.path.dirname(os.path.abspath(__file__))
    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)

    original_dir = os.getcwd()

    if not os.path.isdir(poly_directory):
        os.mkdir(poly_directory)

    os.chdir(poly_directory)

    os.system(this_file_path + "/../polynomial_generation/src/poly-gen_mb-nrg.pl " + str(order) + " " + original_dir + "/" + poly_in_path + " > poly.log")

    os.chdir(original_dir)

def execute_maple(settings, poly_directory):
    """
    Runs maple on the polynomial files to turn them into actual cpp files

    Args:
        settings    - the file containing all relevent settings information
        poly_directory - the directory with all the polynomial files
    """

    this_file_path = os.path.dirname(os.path.abspath(__file__))

    original_dir = os.getcwd()

    os.chdir(poly_directory)

    # clear any old files, because for some reason maple appends to existing files instead of clobbering
    os.system("rm poly-grd.c 2> /dev/null")
    os.system("rm poly-nogrd.c 2> /dev/null")
    os.system("rm poly-grd.cpp 2> /dev/null")
    os.system("rm poly-nogrd.cpp 2> /dev/null")

    os.system("maple poly-grd.maple")
    os.system("maple poly-nogrd.maple")

    os.system(this_file_path + "/../polynomial_generation/src/clean-maple-c.pl < poly-grd.c > poly-grd.cpp")
    os.system(this_file_path + "/../polynomial_generation/src/clean-maple-c.pl < poly-nogrd.c > poly-nogrd.cpp")

    os.chdir(original_dir)

def generate_fit_config(settings, molecule_in, opt_geometry, config_path):
    """
    Generates the fit config for the optimized geometry

    For >=2b, the optimized geometry should be the optimized geometry of each monomer seperated by roughly
    40 angstroms

    Args:
        settings    - the file containing all relevent settings information
        molecule_in     - string of fromat "A1B2"
        opt_geometry - the optimized geometry of this molecule
        config_path - the file to write the config file

    Returns:
        None
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../fitting/src")
    import get_config_data

    if not os.path.isdir(os.path.dirname(config_path)):
        os.mkdir(os.path.dirname(config_path))

    get_config_data.make_config(settings, molecule_in, opt_geometry, config_path)
    
    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../fitting/src")

def generate_1b_fit_code(settings, config, poly_in_path, poly_path, fit_directory):
    """
    Generates the fit code based on the polynomials for a monomer

    Args:
        settings - the file containing all relevent settings information
        config    - monomer config file
        poly_in_path - the A3B2.in type file
        poly_path   - directory where polynomial files are
        fit_directory - directory to generate fit code in

    Returns:
        None
    """

    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../fitting/1B/get_codes")
    import prepare_1b_fitting_code

    if not os.path.isdir(fit_directory):
        os.mkdir(fit_directory)
    
    prepare_1b_fitting_code.prepare_1b_fitting_code(settings, poly_in_path, poly_path, fit_directory)

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../fitting/1B/get_codes") 

def generate_2b_ttm_fit_code(settings, config, molecule_in, fit_directory):
    """
    Generates the fit TTM fit code for a dimer

    Args:
        settings - the file containing all relevent settings information
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
    
    subprocess.call("cp " + this_file_path + "/../fitting/2B/template/* .", shell=True)
    subprocess.call("python " + this_file_path + "/../fitting/2B/get_2b_TTM_codes.py {} {} {}".format(original_dir + "/" + settings, original_dir + "/" + config, molecule_in), shell=True)

    os.chdir(original_dir)   


def compile_fit_code(settings, fit_directory):
    """
    Compiles the fit code in the given directory

    Args:
        settings    - the file containing all relevent settings information
        fit_directory - the directory with the fit code

    Returns:
        None
    """

    original_dir = os.getcwd()

    os.chdir(fit_directory)
    os.system("make clean")
    os.system("make")
    os.chdir(original_dir)

def fit_1b_training_set(settings, fit_code, training_set, fit_directory, fitted_code):
    """
    Fits the fit code to a given training set

    Args:
        settings    - the file containing all relevent settings information
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

    os.system("mv fit-1b.cdl " + fit_directory + "/")
    os.system("mv fit-1b-initial.cdl " + fit_directory + "/")
    os.system("mv correlation.dat " + fit_directory + "/")
    os.system("mv fit-1b.nc " + fitted_code)

def fit_2b_ttm_training_set(settings, fit_code, training_set, fit_directory):
    """
    Fits the ttm fit code to a given training set

    Args:
        settings    - the file containing all relevent settings information
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

        rmds = float(log_lines[-4].split()[2])
        best_rmds = float(best_log_lines[-4].split()[2])

        if rmds < best_rmds:
            os.system("mv " + fit_directory + "/fit.log " + fit_directory + "/best_fit.log")
            

        attempts += 1

    os.system("rm " + fit_directory + "/fit.log")
    os.system("mv individual_terms.dat " + fit_directory + "/")
    os.system("mv ttm-params.txt " + fit_directory + "/")
    os.system("mv correlation.dat " + fit_directory + "/")
