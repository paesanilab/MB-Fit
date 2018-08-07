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
        'files':            {'log_path': 'outputs'},
        'molecule':         {'charges': '0',
                             'spins': '1'},
        'config_generator': {'optimize': 'true',
                             'random': 'P',
                             'geometric': 'false',
                             'linear': 'true', 
                             'code': 'psi4'},
        'psi4':             {'memory': '1GB',
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

    output_id = name + "_" + config['config_generator']['code'] + "_" + config['config_generator']['method'] + "_" + config['config_generator']['basis']
    log_name = config['files']['log_path'] + "/" + output_id

    filenames = {
        "input_geometry": input_geo_fn,
        "optimized_geometry": log_name + "_optimized.xyz",
        "normal_modes": log_name + "_normalmodes.dat",
        "gcn_input": log_name + "_gcn.inp",
        "gcn_output": log_name + "_gcn.out",
        "norm_config": log_name + "_configurations.xyz",
        "qchem_opt_input": log_name + "_qchem_optimization.inp",
        "qchem_freq_input": log_name + "_qchem_frequencies.inp",
        "qchem_opt_output": log_name + "_qchem_optimization.out",
        "qchem_freq_output": log_name + "_qchem_frequencies.out"
    }

    if "optimized_geometry" in config["files"]:
        filenames["optimized_geometry"] = config["files"]["optimized_geometry"]

    if "xyz_files" in config["files"]:
        filenames["norm_config"] = config["files"]["xyz_files"] + "/" + output_id + "_configs.xyz"
    
    return filenames, log_name
