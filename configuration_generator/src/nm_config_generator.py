# normalconfigurationgenerator
# 
# Wrapper that inputs values into Sandra's FORTRAN scripts.
#
# @author Ronak

import os

def generate(molecule, gcn_input_path, gcn_output_path, optimized_geometry_path, normal_modes_path, norm_config_path, gcn_executor_path, dim_null, random, num_configs, geometric, linear):
    print("Running normal distribution configuration generator...")
    with open(gcn_input_path, 'w') as input_file:
        input_file.write("\'" + optimized_geometry_path + "\'\n")             # optimized geometry
        input_file.write("\'" + normal_modes_path + "\'\n")                 # normal modes
        input_file.write(str(3 * molecule.natom()) + " " + str(dim_null) + "\n")    # dim; dimnull
        input_file.write(random + " " + str(num_configs) + "\n")            # random method; nconfigs
        input_file.write("\'" + norm_config_path + "\'\n")                    # output one-body configurations
        input_file.write(geometric + " " + linear + "\n")                    # geometric, linear
        input_file.write(".true.")                                            # verbose
        
#    os.system("bash gcn/src/generate_configs_normdistrbn < " + gcn_input_path + " > " + gcn_output_path)
    
    with open(gcn_executor_path, 'w') as executor_file:
        executor_file.write("#!/bin/bash\n\n")
        executor_file.write("generate_configs_normdistrbn < " + gcn_input_path + " > " + gcn_output_path + "\n")
        executor_file.write("exit 0")
        
    os.system("bash " + gcn_executor_path)
    
    print("Normal Distribution Configuration generation complete.")

