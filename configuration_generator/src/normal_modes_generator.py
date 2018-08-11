import configparser
import qcalc
import output_writer
import molecule
import psi4_helper

def generate_normal_modes(settings, geo_path, normal_modes_file):
    
    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)

    psi4_helper.init(config, config["files"]["log_path"] + "/" + "optimize")

    with open(geo_path, "r") as geo_file:
        molec = molecule.Molecule(geo_file.read())
    
    normal_modes, frequencies, red_masses = qcalc.frequencies(molec, {"log_name": config["files"]["log_path"]}, config)
    dim_null = 3 * molec.num_atoms - len(normal_modes)
    
    output_writer.write_normal_modes(normal_modes, frequencies, red_masses, normal_modes_file)
    print("DIM NULL: " + str(dim_null))
