# Psi4Helper.py
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

import psi4helper
import config_loader
import nm_config_generator

if len(sys.argv) == 2:
    parse_fail, name, input_geo, linear_geo, charge, multiplicity, random, num_configs, geometric, linear, method, basis = config_loader.parse_config(sys.argv[1])
else:
    parse_fail = True
    print("Usage: automated-psi4.py [config.txt]")
    
if parse_fail:
    print("ERROR: The configuration file does not conform to the input format.")
else:
    # filepaths
    optimized_geometry_path = "inputs/" + name + "_psi4_" + method + "_" + basis + "_optimized.xyz"
    normal_modes_path = "inputs/" + name + "_psi4_" + method + "_" + basis + "_normalmodes.dat"

    gcn_input_path = "gcn/" + name + "_psi4_" + method + "_" + basis + ".inp"
    gcn_output_path = "gcn/" + name + "_psi4_" + method + "_" + basis + ".out"

    gcn_executor_path = "gcn/" + name + "_psi4_" + method + "_" + basis + ".sh"

    norm_config_path = "gcn/" + name + "_psi4_" + method + "_" + basis + "_configurations.xyz"

    # psi4 Global Parameters
    memory = "1GB"
    num_threads = 2
    output_path = "logs/" + name + "_psi4_" + method + "_" + basis + ".out"


    psi4.core.set_output_file(output_path, False)
    psi4.set_memory(memory)
    psi4.set_num_threads(num_threads)

    # Defining the molecule
    geo = input_geo + "\n" + charge + " " + multiplicity
    
    mol = psi4.geometry(geo)

    # Defining the computation models for steps 1 and 2
    model = method + "/" + basis
    
    if linear_geo:
        dim_null = 5
    else:
        dim_null = 6
    
    # Step 1
    psi4helper.optimize_geometry(mol, model, optimized_geometry_path)
    
    # Step 2
    psi4helper.frequency_calculations(mol, model, dim_null, normal_modes_path)
    
    # Step 3
    nm_config_generator.generate(mol, gcn_input_path, gcn_output_path, optimized_geometry_path, normal_modes_path, norm_config_path, gcn_executor_path, dim_null, random, num_configs, geometric, linear)

