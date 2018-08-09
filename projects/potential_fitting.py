# main driver for potential fitting program, takes only a directory name

import sys
import os
import configparser

def potential_fit(project_directory):
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

    if not os.path.isdir(config['files']['log_path']):
        os.mkdir(config['files']['log_path'])
    if not os.path.isdir(config['files']['xyz_files']):
        os.mkdir(config['files']['xyz_files'])
    if not os.path.isdir(config['files']['poly_path']):
        os.mkdir(config['files']['poly_path'])
    if not os.path.isdir(config['files']['training_set']):
        os.mkdir(config['files']['training_set'])

    # check spin and set energy calculation method (defaulted to HF)
    
    if not config['molecule']['spins'] == 1:
        config.set('energy_calculator', 'method', 'UHF')
        print ("changed method to " + config['energy_calculator']['method'])

    # Step 1: generate configurations

    generate_configs(project_directory, config)

    # Step 2: initialize database

    initialize_database(project_directory, config)

    # Step 3: fill database

    fill_database(project_directory, config)

    # Step 4: create input file for polynomial

    generate_input(project_directory, config)

    # Step 5: generate polynomials

    os.chdir(config['files']['poly_path'])

    generate_polynomials(project_directory, config)

    # Step 6: run maple on polynomails

    #run_maple(project_directory, config)

    os.chdir("..")

    # Step 7: generate training set xyz file

    generate_training_set(project_directory, config)

    # Step 8: generate fitting code

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

def initialize_database(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../qm_mb_energy_calculator/src/database_initializer.py settings.ini " + config['files']['database'] + " " + config['files']['xyz_files'])

def fill_database(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../qm_mb_energy_calculator/src/database_filler.py settings.ini " + config['files']['database'] + " " + config['files']['xyz_files'])

def generate_input(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../polynomial_generation/src/generate_input_poly.py " + config['files']['poly_in_path'])

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
      

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Incorrect number of command line arguments")
        print("Usage: python potential_fitting.py <project_directory>")
        exit(1)
    potential_fit(sys.argv[1])
