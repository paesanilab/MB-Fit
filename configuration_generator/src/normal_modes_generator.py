import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../qm_mb_energy_calculator/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import qcalc
import molecule_parser
import settings_reader

def generate_normal_modes(settings_path, geo_path, normal_modes_file):
    
    settings = settings_reader.SettingsReader(settings_path)

    with open(geo_path, "r") as geo_file:
        molecule = molecule_parser.xyz_to_molecules(geo_path, settings)[0]
    
    normal_modes, frequencies, red_masses = qcalc.frequencies(settings, molecule, settings.get("config_generator", "method"), settings.get("config_generator", "basis"))

    dim_null = 3 * molecule.get_num_atoms() - len(normal_modes)
    
    write_normal_modes(normal_modes, frequencies, red_masses, normal_modes_file)

    return dim_null

def write_normal_modes(normal_modes, frequencies, red_masses, normal_modes_fn):

    norm_formatter = "{:> 12.4f}"
    freq_formatter = "{:.2f}"
    mass_formatter = "{:.6f}"

    normal_out = ""
    
    num_modes = len(frequencies)

    for i in range(1, 1 + num_modes):
        index = i - 1
        
        geo = normal_modes[index]
        
        normal_out += "normal mode: " + str(i) + "\n"
        normal_out += freq_formatter.format(frequencies[index]) + "\n"
        normal_out += "red_mass = " + mass_formatter.format(red_masses[index]) + "\n"
        
        for atom in range(len(geo)):
            normal_out += norm_formatter.format(0.0 if abs(float(geo[atom][0])) < 1e-6 else float(geo[atom][0])) + "\t"
            normal_out += norm_formatter.format(0.0 if abs(float(geo[atom][1])) < 1e-6 else float(geo[atom][1])) + "\t"
            normal_out += norm_formatter.format(0.0 if abs(float(geo[atom][2])) < 1e-6 else float(geo[atom][2])) + "\n"
        
        normal_out += "\n"

    with open(normal_modes_fn, 'w') as norm_file:
        norm_file.write(normal_out)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python geometry_optimizer.py <settings> <geo_path> <normal_modes_path>")
        exit(1)
    generate_normal_modes(sys.argv[1], sys.argv[2], sys.argv[3])
