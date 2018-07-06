# gcn_runner.py
# 
# Wrapper that calls Sandra's FORTRAN GCN executable.
#
# @author Ronak

import os

def generate(config, dim_null, num_atoms, filenames):
    print("Running normal distribution configuration generator...")

    gcn_params = config['config_generator']

    with open(filenames['gcn_input'], 'w') as input_file:
        input_file.write("\'" + filenames['optimized_geometry'] + "\'\n") # optimized geometry
        input_file.write("\'" + filenames['normal_modes'] + "\'\n") # normal modes
        input_file.write(str(3 * num_atoms) + " " + str(dim_null) + "\n") # dim; dimnull
        input_file.write(gcn_params['random'] + " " + gcn_params['num_configs'] + "\n") # random method; nconfigs
        input_file.write("\'" + filenames['norm_config'] + "\'\n") # output one-body configurations
        input_file.write(gcn_params['geometric'] + " " + gcn_params['linear'] + "\n") # geometric, linear
        input_file.write(".true.") # verbose
        
    
    script_path = os.path.dirname(os.path.realpath(__file__))

    path = script_path[:-4]
    os.system(path + "/norm_distribution/src/generate_configs_normdistrbn < " + filenames['gcn_input'] + " > " + filenames['gcn_output'])
    
    print("Normal Distribution Configuration generation complete.")

