import configparser
import molecule
import os

def generate_configurations(settings, geo_path, normal_modes, dim_null, config_dir):
    print("Running normal distribution configuration generator...")
    
    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)
    
    with open(geo_path, "r") as geo_file:
        molec = molecule.Molecule(geo_file.read())


    with open(config["files"]["log_path"] + "/" + "config_input.log", "w") as input_file:
        input_file.write("'" + geo_path + "'\n") # optimized geometry
        input_file.write("'" + normal_modes + "'\n") # normal modes
        input_file.write(str(3 * molec.num_atoms) + " " + str(dim_null) + "\n") # dim; dimnull
        input_file.write(config["config_generator"]["random"] + " " + config["config_generator"]["num_configs"] + "\n") # random method; nconfigs
        input_file.write("'" + config_dir + "/configs.xyz'\n") # output one-body configurations
        input_file.write(config["config_generator"]["geometric"] + " " + config["config_generator"]["linear"] + "\n") # geometric, linear
        input_file.write(".true.") # verbose
    
    os.system(os.path.dirname(os.path.abspath( __file__ )) + "/../norm_distribution/src/generate_configs_normdistrbn < " + config["files"]["log_path"] + "/config_input.log" + " > " + config["files"]["log_path"] + "/config_output.log")

    os.system("cp " + geo_path + " " + config_dir + "/geo.opt.xyz")
    
    print("Normal Distribution Configuration generation complete.")
        
