# main driver for potential fitting program, takes only a directory name

import sys
import os
import configparser

def potential_fit(project_directory, start_location = 0):
    """
    Runs the whole potential fitting workflow on a project directory
    """

    # first make sure a settings.ini exists in project_directory
    config = configparser.ConfigParser(allow_no_value=False)

    try:
        config.read(project_directory + "/settings.ini")
        config.set("files", "directory", os.path.dirname(os.path.abspath(__file__)))
    except:
        print("settings.ini file in project directory either missing or unopenable")
        exit(1)

    # move into directory
    os.chdir(project_directory)

    # create all required directories

    if start_location <= 0:
        if not os.path.isdir(config['files']['log_path']):
            os.mkdir(config['files']['log_path'])
        if not os.path.isdir(config['files']['xyz_files']):
            os.mkdir(config['files']['xyz_files'])
        if not os.path.isdir(config['files']['poly_path']):
            os.mkdir(config['files']['poly_path'])
        if not os.path.isdir(config['files']['training_set']):
            os.mkdir(config['files']['training_set'])

    # Step 1: generate configurations

    if start_location <= 1:
        generate_configs(project_directory, config)

    # Step 2: initialize database

    if start_location <= 2:
        initialize_database(project_directory, config)

    # Step 3: fill database

    if start_location <= 3:
        fill_database(project_directory, config)

    # Step 4: create input file for polynomial

    if start_location <= 4:
        generate_input(project_directory, config)

    # Step 5: generate polynomials

    if start_location <= 5:
        os.chdir(config['files']['poly_path'])
        generate_polynomials(project_directory, config)
        os.chdir("..")

    # Step 6: generate training set xyz file

    if start_location <= 6:
        generate_training_set(project_directory, config)

    # Step 7: generate fitting code

    #os.chdir(config['files']['fit_code'])

    # generate_fitting(project_directory, config)

    #os.chdir("..")

    # Step 9: fit the code

    # remove unneeded empty directories
    os.system("find . -type d -empty -delete")

"""
YES: the system calls are kinda messy code. But at the moment, it is the best way to do it, because many of these modules do not have methods that can be called.

while the database stuffs can be called as a method rather than from the command line, the configuration generator (at least as I am writing this) must be
called from the command line. So I called everything from the command line for consistancies sake
"""

def generate_configs(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../configuration_generator/src/nmcgen.py settings.ini")
    os.system("cp " + config['files']['optimized_geometry'] + " " + config['files']['xyz_files'] + "/config.opt.xyz")

def initialize_database(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../qm_mb_energy_calculator/src/database_initializer.py settings.ini " + config['files']['database'] + " " + config['files']['xyz_files'])

def fill_database(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../qm_mb_energy_calculator/src/database_filler.py settings.ini " + config['files']['database'] + " " + config['files']['xyz_files'])

def generate_input(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../polynomial_generation/src/generate_input_poly.py settings.ini " + config['files']['poly_in_path'])

def generate_polynomials(project_directory, config):
    os.system(config['files']['directory'] + "/../polynomial_generation/src/poly-gen_mb-nrg.pl " + config['poly_generator']['order'] + " " + config['files']['directory'] + "/" + project_directory + "/" + config['files']['poly_in_path'] + " > " + config['files']['directory'] + "/" + project_directory + "/" + config['files']['poly_path'] + "/poly.log")

def run_maple(project_directory, config):
    # maple commands can only run on maple machine

    os.system("maple poly-grd.maple")
    os.system("maple poly-nogrd.maple")

    os.system(config['files']['directory'] + "/../polynomial_generation/src/clean-maple-c.pl < poly-grd.c > poly-grd.cpp")
    os.system(config['files']['directory'] + "/../polynomial_generation/src/clean-maple-c.pl < poly-nogrd.c > poly-nogrd.cpp")

def generate_training_set(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../qm_mb_energy_calculator/src/database_reader.py settings.ini " + config['files']['database'] + " " + config['files']['training_set'] + "/training_set.xyz")
    
def generate_fitting(project_directory, config):

    # copy over A1B2.in file
    os.system("cp " + config['files']['directory'] + "/" + project_directory + "/" + config['files']['poly_in_path'] + " .")
    
    # copy over the poly-direct.cpp file
    os.system("cp " + config['files']['directory'] + "/" + project_directory + "/" + config['files']['poly_path'] + "/poly-direct.cpp .")

    # copy over and rename the poly-grd.cpp file
    os.system("cp " + config['files']['directory'] + "/" + project_directory + "/" + config['files']['poly_path'] + "/poly-grd.cpp ./poly_1b_" + config['files']['poly_in_path'][:-3] + "_v1x.cpp")

    # copy over and rename the poly-nogrd.cpp file
    os.system("cp " + config['files']['directory'] + "/" + project_directory + "/" + config['files']['poly_path'] + "/poly-nogrd.cpp ./poly_1b_" + config['files']['poly_in_path'][:-3] + "_v1.cpp")

    # copy over and rename poly-model.h file
    os.system("cp " + config['files']['directory'] + "/" + project_directory + "/" + config['files']['poly_path'] + "/poly-model.ch ./poly_1b_" + config['files']['poly_in_path'][:-3] + "_v1x.h")

    
    
    os.system("python " + config['files']['directory'] + "/../fitting/1B/get_codes/prepare_1b_fitting_code.sh" + config['files']['directory'] + "/" + project_directory + config['files']['poly_in_path'] + " " + config['files']['directory'] + "/" + config['files']['poly_path'])
      

import sys
import os
import configparser

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
        None

    Other Effects:
        Writes DIM NULL to console
    """
    
    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../configuration_generator/src")
    import normal_modes_generator
   
    if not os.path.isdir(os.path.dirname(normal_modes)):
        os.mkdir(os.path.dirname(normal_modes))

    normal_modes_generator.generate_normal_modes(settings, geo, normal_modes)

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../configuration_generator/src") 


def generate_configurations(settings, geo, normal_modes, dim_null, configurations):
    """
    Generates configurations for a given molecule from a set of normal modes

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
    print(sys.path)
    import database_initializer

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

    database_filler.fill_database(settings, database_name, "unused")

    sys.path.remove(os.path.dirname(os.path.abspath(__file__)) + "/../qm_mb_energy_calculator/src") 

def generate_training_set(settings, database_name, training_set, method = "%", basis = "%", cp = "%"):
    """
    Generates a training set from the energies inside a database.

    Specific method, bnasis, and cp may be specified to only use energies calculated
    with a specific model.

    '%' can be used to stand in as a wild card, meaning any method/basis/cp is ok

    Args:
        settings    - the file containing all relevent settings information
        database_name - the file in which the database is stored
        training_set - the file to write the training set to
        method      - only use energies calcualted by this method
        basis       - only use energies calculated in this basis
        cp          - only use energies calculated with the same cp

    Returns:
        None
    """

    # imports have to be here because python is bad
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../qm_mb_energy_calculator/src")
    import training_set_generator

    if not os.path.isdir(os.path.dirname(training_set)):
        os.mkdir(os.path.dirname(training_set))

    training_set_generator.generate_training_set(settings, database_name, training_set, method, basis, cp)

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

    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)

    original_dir = os.getcwd()

    if not os.path.isdir(poly_directory):
        os.mkdir(poly_directory)

    os.chdir(poly_directory)

    os.system(os.path.dirname(os.path.abspath(__file__)) + "/../polynomial_generation/src/poly-gen_mb-nrg.pl " + str(order) + " " + original_dir + "/" + poly_in_path + " > poly.log")

    os.chdir(original_dir)

def execute_maple(settings, poly_directory):
    """
    Runs maple on the polynomial files to turn them into actual cpp files

    Args:
        settings    - the file containing all relevent settings information
        poly_directory - the directory with all the polynomial files
    """

    original_dir = os.getcwd()

    os.chdir(poly_directory)

    os.system("maple poly-grd.maple")
    os.system("maple poly-nogrd.maple")

    os.system(os.path.dirname(os.path.abspath(__file__)) + "/../polynomial_generation/src/clean-maple-c.pl < poly-grd.c > poly-grd.cpp")
    os.system(os.path.dirname(os.path.abspath(__file__)) + "/../polynomial_generation/src/clean-maple-c.pl < poly-nogrd.c > poly-nogrd.cpp")

    os.chdir(original_dir)

def generate_fit_code(settings, poly_in_path, poly_path, fit_directory):
    """
    Generates the fit code based on the polynomials

    Only works for 1b right now

    Args:
        settings    - the file containing all relevent settings information
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

def fit_training_set(settings, fit_code, training_set, fitted_code):
    """
    Fits the fit code to a given training set

    Args:
        settings    - the file containing all relevent settings information
        fit_code    - the code to fit
        training_set - the training set to fit the code to
        fitted_code - file to write final fitted code to

    Returns:
        None
    """
    
    os.system(fit_code + " " + training_set)

if __name__ == "__main__":
    if len(sys.argv) == 2:
        potential_fit(sys.argv[1])
    elif len(sys.argv) == 3:
        potential_fit(sys.argv[1], int(sys.argv[2]))
    else:
        print("Incorrect number of command line arguments")
        print("Usage: python potential_fitting.py <project_directory>")
        exit(1)
