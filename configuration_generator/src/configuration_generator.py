import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../qm_mb_energy_calculator/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import molecule_parser
import settings_reader
def generate_configurations(settings_path, geo_path, normal_modes, dim_null, config_path):
    print("Running normal distribution configuration generator...")
    
    settings = settings_reader.SettingsReader(settings_path)
    
    with open(geo_path, "r") as geo_file:
        molecule = molecule_parser.xyz_to_molecules(geo_path, settings)[0]

    input_path = settings.get("files", "log_path") + "/config_generator/{}.in".format(molecule.get_name())

    if not os.path.exists(os.path.dirname(input_path)):
        os.makedirs(os.path.dirname(input_path))

    with open(input_path, "w") as input_file:
        input_file.write("'" + geo_path + "'\n") # optimized geometry
        input_file.write("'" + normal_modes + "'\n") # normal modes
        input_file.write(str(3 * molecule.get_num_atoms()) + " " + str(dim_null) + "\n") # dim; dimnull
        input_file.write(settings.get("config_generator", "random") + " " + settings.get("config_generator", "num_configs") + "\n") # random method; nconfigs
        input_file.write("'" + config_path + "'\n") # output one-body configurations
        input_file.write(settings.get("config_generator", "geometric") + " " + settings.get("config_generator", "linear") + "\n") # geometric, linear
        input_file.write(".true.") # verbose
    
    log_path = settings.get("files", "log_path") + "/config_generator/{}.log".format(molecule.get_name())
    

    os.system(os.path.dirname(os.path.abspath( __file__ )) + "/../norm_distribution/src/generate_configs_normdistrbn < " + input_path + " > " + log_path)

    print("Normal Distribution Configuration generation complete.")
        
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python geometry_optimizer.py <settings> <geo_path> <normal_modes_path> <dim null> <config_path>")
        exit(1)
    generate_configurations(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5])
