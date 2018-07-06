# config_loader.py
#
# Loads the configuration .ini file for nmcgen

from configparser import ConfigParser
from argparse import ArgumentParser

def load():
    arg_parser = ArgumentParser(description='Generates configurations based on normal modes.')
    arg_parser.add_argument('config', metavar='config', type=str, nargs=1, help='Provide the path to the .ini containing the settings for this script.')

    args = arg_parser.parse_args()

    config = ConfigParser(allow_no_value=False)

    # Config Default Values
    config.read_dict({
        'files':            {'output_path': 'outputs/'},
        'molecule':         {'charge': '0',
                             'multiplicity': '1'},
        'config_generator': {'random': 'P',
                             'geometric': 'false',
                             'linear': 'true'},
        'program':          {'code': 'psi4',
                             'memory': '1GB',
                             'num_threads': '1'}
    })

    config.read(args.config[0])
    
    return config
    
def process_files(config):
    # filenames
    input_geo_fn = config['files']['input_geometry']

    if input_geo_fn.endswith('.xyz'):
        slash_split = input_geo_fn.split("/")
        
        name = slash_split[len(slash_split) - 1]
        name = input_geo_fn[:-4]
        name = name.split("_")[0]
    else:
        raise IOError("Input geometry not an .xyz")  

    if 'name' in config['files']:
        name = config['files']['name']

    output_id = name + "_" + config['program']['code'] + "_" + config['model']['method'] + "_" + config['model']['basis']
    log_name = config['files']['output_path'] + output_id

    filenames = {
        "optimized_geometry": log_name + "_optimized.xyz",
        "normal_modes": log_name + "_normalmodes.dat",
        "gcn_input": log_name + "_gcn.inp",
        "gcn_output": log_name + "_gcn.out",
        "norm_config": log_name + "_configurations.xyz",
    }
    
    return filenames, log_name
