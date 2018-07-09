# output_writer.py
#
# Functions for writing optimized geometries or normal modes for gcn


opt_formatter =  "{:< 12.6f}"
energy_formatter = "{:.6f}"

norm_formatter = "{:> 12.4f}"
freq_formatter = "{:.2f}"
mass_formatter = "{:.6f}"

def write_optimized_geo(molecule, energy, optimized_geometry_fn):
    with open(optimized_geometry_fn, 'w') as opt_geo_file:
            opt_geo_file.write(str(molecule.num_atoms) + "\n")
            opt_geo_file.write(energy_formatter.format(energy) + "\n")
            
            for i in range(molecule.num_atoms):
                atom = molecule.atoms[i]
                coordinates = molecule.coordinates[i]
                
                x = filter_neg_zero(opt_formatter.format(coordinates[0]), 6) 
                y = filter_neg_zero(opt_formatter.format(coordinates[1]), 6)
                z = filter_neg_zero(opt_formatter.format(coordinates[2]), 6)
                
                opt_geo_file.write(atom + "\t" + x + " " + y + " " + z + "\n")
                
def write_normal_modes(normal_modes, frequencies, red_masses, normal_modes_fn):
    normal_out = ""
    
    num_modes = len(frequencies)

    for i in range(1, 1 + num_modes):
        index = i - 1
        
        geo = normal_modes[index]
        
        normal_out += "normal mode: " + str(i) + "\n"
        normal_out += freq_formatter.format(frequencies[index]) + "\n"
        normal_out += "red_mass = " + mass_formatter.format(red_masses[index]) + "\n"
        
        for atom in range(len(geo)):
            normal_out += filter_neg_zero(norm_formatter.format(float(geo[atom][0])), 4) + "\t"
            normal_out += filter_neg_zero(norm_formatter.format(float(geo[atom][1])), 4) + "\t"
            normal_out += filter_neg_zero(norm_formatter.format(float(geo[atom][2])), 4) + "\n"
        
        normal_out += "\n"

    with open(normal_modes_fn, 'w') as norm_file:
        norm_file.write(normal_out)
                
def filter_neg_zero(string, num_zeros):
    zeros = num_zeros * "0"
    
    if "-0." + zeros in string:
        string = str.replace(string, "-0." + zeros, " 0." + zeros)
        
    return string
