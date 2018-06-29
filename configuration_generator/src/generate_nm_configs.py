# generate_nm_configs.py
# 
# Automation of the following steps of the workflow:
#	1. geometry optimization (using psi4)
#	2. generation of normal modes/frequency calculation (using psi4)
#	3. configuration generation (using Sandra's code)
# Based on the user input of the initial geometry and the calculation model
#
# @author Ronak

import psi4
import os
import sys

from configparser import ConfigParser
from argparse import ArgumentParser

import psi4_helper
import nm_config_generator

arg_parser = ArgumentParser(description='Generates the normal configurations.')
arg_parser.add_argument('config', metavar='config', type=str, nargs=1, help='Provide the path the the .ini containing the parameters to run the script.')

args = arg_parser.parse_args()

config = ConfigParser(allow_no_value=False)
config.read(args.config[0])

# filepaths
input_geo_path = config['file_input']['input_geometry_path']

if input_geo_path.endswith('.xyz'):
    name = input_geo_path[:-4]
else:
    raise IOError("Input geometry not an .xyz")

input_geo = ""
with open(input_geo_path, 'r') as input_file:
    for i, line in enumerate(input_file):
        if i > 1:
            input_geo += line    

log_name = name + "_" + config['program']['code'] + "_" + config['model']['method'] + "_" + config['model']['basis']

optimized_geometry_path = "outputs/" + log_name + "_optimized.xyz"
normal_modes_path = "outputs/" + log_name + "_normalmodes.dat"
    
gcn_input_path = "outputs/" + log_name + "_gcn.inp"
gcn_output_path = "outputs/" + log_name + "_gcn.out"
    
norm_config_path = "outputs/" + log_name + "_configurations.xyz"

output_path = "outputs/" + log_name + ".log"

if config['program']['code'] == "psi4":
    psi4.core.set_output_file(output_path, False)
    psi4.set_memory(config['program']['memory'])
    psi4.set_num_threads(int(config['program']['num_threads']))

    # Defining the molecule
    geo = input_geo + "\n" + config['molecule']['charge'] + " " + config['molecule']['multiplicity']
    mol = psi4.geometry(geo)

    # Step 1
    psi4_helper.optimize_geometry(mol, config, optimized_geometry_path)

    # Step 2
    dim_null = psi4_helper.frequency_calculations(mol, config, normal_modes_path)

    # Step 3
    nm_config_generator.generate(mol, gcn_input_path, gcn_output_path, optimized_geometry_path, normal_modes_path, norm_config_path, config, dim_null)

