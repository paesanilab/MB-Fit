# nm_config_generator.py
# 
# Wrapper that inputs values into Sandra's FORTRAN scripts.
#
# @author Ronak

import os

def generate(molecule, gcn_input_path, gcn_output_path, optimized_geometry_path, normal_modes_path, norm_config_path, gcn_executor_path, config, dim_null):
    print("Running normal distribution configuration generator...")

    cg_params = config['config_generator']

    with open(gcn_input_path, 'w') as input_file:
        input_file.write("\'" + optimized_geometry_path + "\'\n")					# optimized geometry
        input_file.write("\'" + normal_modes_path + "\'\n")							# normal modes
        input_file.write(str(3 * molecule.natom()) + " " + str(dim_null) + "\n")    # dim; dimnull
        input_file.write(cg_params['random'] + " " + str(3 * molecule.natom() - dim_null) + "\n")       # random method; nconfigs
        input_file.write("\'" + norm_config_path + "\'\n")                    		# output one-body configurations
        input_file.write(cg_params['geometric'] + " " + cg_params['linear'] + "\n") # geometric, linear
        input_file.write(".true.")                                            		# verbose
        
    os.system("./../norm_distribution/src/generate_configs_normdistrbn < " + gcn_input_path + " > " + gcn_output_path)
    
#    with open(gcn_executor_path, 'w') as executor_file:
#        executor_file.write("#!/bin/bash\n\n")
#        executor_file.write("../norm_distribution/src/generate_configs_normdistrbn < " + gcn_input_path + " > " + gcn_output_path + "\n")
#        executor_file.write("exit 0")
        
#    os.system("bash " + gcn_executor_path)
    
    print("Normal Distribution Configuration generation complete.")

