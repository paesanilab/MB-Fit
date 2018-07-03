# main driver for potential fitting program, takes only a directory name

import sys
import os
import configparser

def potential_fit(project_directory):

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

    # create the logs directory for log files

    if not os.path.isdir(config['files']['log_path']):
        os.mkdir(config['files']['log_path'])
    if not os.path.isdir(config['files']['xyz_files']):
        os.mkdir(config['files']['xyz_files'])
    if not os.path.isdir(config['files']['poly_path']):
        os.mkdir(config['files']['poly_path'])
    if not os.path.isdir(config['files']['training_set']):
        os.mkdir(config['files']['training_set'])

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

    os.chdir("..")

    # Step 6: generate training set xyz file

    generate_training_set(project_directory, config)

    os.system("find . -type d -empty -delete")

def generate_configs(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../configuration_generator/src/nmcgen.py settings.ini")

def initialize_database(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../qm_mb_energy_calculator/src/database_initializer.py settings.ini " + config['files']['database'] + " " + config['files']['xyz_files'])

def fill_database(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../qm_mb_energy_calculator/src/database_filler.py settings.ini " + config['files']['database'] + " " + config['files']['xyz_files'])

def generate_input(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../polynomial_generation/src/generate_input_poly.py " + config['files']['poly_in_path'])

def generate_polynomials(project_directory, config):
    os.system(config['files']['directory'] + "/../polynomial_generation/src/poly-gen_mb-nrg.pl " + config['poly_generator']['order'] + " " + config['files']['directory'] + "/" + project_directory + "/" + config['files']['poly_in_path'] + " > " + config['files']['directory'] + "/" + project_directory + "/" + config['files']['log_path'] + "/poly_log.log")  

def generate_training_set(project_directory, config):
    os.system("python " + config['files']['directory'] + "/../qm_mb_energy_calculator/src/database_reader.py settings.ini " + config['files']['database'] + " " + config['files']['training_set'] + "/training_set.xyz")
    
      

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Incorrect number of command line arguments");
        print("Usage: python potential_fitting.py <project_directory>");
        exit(1)
    potential_fit(sys.argv[1])
