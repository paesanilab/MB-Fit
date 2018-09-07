import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../qm_mb_energy_calculator/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import qcalc
import output_writer
import molecule_parser
import settings_reader

def generate_normal_modes(settings_path, geo_path, normal_modes_file):
    
    settings = settings_reader.SettingsReader(settings_path)

    with open(geo_path, "r") as geo_file:
        molecule = molecule_parser.xyz_to_molecules(geo_path, settings)[0]
    
    normal_modes, frequencies, red_masses = qcalc.frequencies(settings, molecule, settings.get("config_generator", "method"), settings.get("config_generator", "basis"))

    dim_null = 3 * molecule.get_num_atoms() - len(normal_modes)
    
    output_writer.write_normal_modes(normal_modes, frequencies, red_masses, normal_modes_file)

    return dim_null

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python geometry_optimizer.py <settings> <geo_path> <normal_modes_path>")
        exit(1)
    generate_normal_modes(sys.argv[1], sys.argv[2], sys.argv[3])
