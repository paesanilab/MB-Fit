# nmcgen.py
# 
# Normal Mode Configuration Generator
#
# Automation of the following steps of the workflow:
#	1. geometry optimization
#	2. generation of normal modes/frequency calculation
#	3. configuration generation
#
# The calculations for steps one and two are completed using a quantum chemistry code, psi4 by default. Step three is completed using Sandra's GCN code.
#
# @author Ronak

import psi4
import os
import sys

from configparser import ConfigParser
from argparse import ArgumentParser

import psi4_helper
import gcn_runner

arg_parser = ArgumentParser(description='Generates configurations based on normal modes.')
arg_parser.add_argument('config', metavar='config', type=str, nargs=1, help='Provide the path to the .ini containing the settings for this script.')

args = arg_parser.parse_args()

config = ConfigParser(allow_no_value=False)

# Config Default Values
config.read_dict({
    'files':            {'output_path': 'outputs/'},
    'molecule':         {'charge': '0',
                         'multiplicity': '1'},
    'config_generator': {'random': 'Q',
                         'geometric': 'false',
                         'linear': 'true'},
    'program':          {'code': 'psi4',
                         'memory': '1GB',
                         'num_threads': '1'}
})

config.read(args.config[0])

# filenames
input_geo_fn = config['files']['input_geometry']

if input_geo_fn.endswith('.xyz'):
    name = input_geo_fn[:-4]
else:
    raise IOError("Input geometry not an .xyz")

input_geo = ""
with open(input_geo_fn, 'r') as input_file:
    for i, line in enumerate(input_file):
        if i > 1:
            input_geo += line    

output_id = name + "_" + config['program']['code'] + "_" + config['model']['method'] + "_" + config['model']['basis']
log_name = config['files']['output_path'] + output_id

filenames = {"optimized_geometry": log_name + "_optimized.xyz",
             "normal_modes": log_name + "_normalmodes.dat",
             "gcn_input": log_name + "_gcn.inp",
             "gcn_output": log_name + "_gcn.out",
             "norm_config": log_name + "_configurations.xyz",
}

if config['program']['code'] == "psi4":
    psi4.core.set_output_file(log_name + ".log", False)
    psi4.set_memory(config['program']['memory'])
    psi4.set_num_threads(int(config['program']['num_threads']))

    # Defining the molecule
    geo = input_geo + "\n" + config['molecule']['charge'] + " " + config['molecule']['multiplicity']
    mol = psi4.geometry(geo)

    # Step 1
    psi4_helper.optimize_geometry(mol, config, filenames['optimized_geometry'])

    # Step 2
    dim_null = psi4_helper.frequency_calculations(mol, config, filenames['normal_modes'])

    # Step 3
    gcn_runner.generate(mol, config, dim_null, filenames)

