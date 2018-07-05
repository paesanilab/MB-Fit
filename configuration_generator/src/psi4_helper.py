# psi4_helper
# 
# Allows for function calls to the Psi4 library in order to optimize the molecular geometry and to generate the normal modes/caluclate their frequencies.
#
# @author Ronak

import psi4
import os
import sys

# Helper functions

# Converts an imaginary output from psi4 to a real output by splicing off the +j term (Sandra's code currently doesn't input imaginary numbers, so only the real frequency portion is maintained). I'm not sure if this is chemically/physically/mathematically correct, but it's what I'm doing.
def treat_imaginary(contents):
    string = str(contents)
    
    result = string.replace("(","")
    result = result.replace(")","")
    
    plus_loc = result.index("+")
    result = result[:plus_loc]
    
    return result
    
def filter_neg_zero(string, num_zeros):
    zeros = num_zeros * "0"
    
    if "-0." + zeros in string:
        string = str.replace(string, "-0." + zeros, " 0." + zeros)
        
    return string

# Main Functions
def optimize(molecule, config, optimized_geometry_fn):
    print("Optimizing geometry...")

    opt_formatter =  "{:< 12.6f}"
    energy_formatter = "{:.6f}"

    e = psi4.optimize(config['config_generator']['method'] + '/' + config['config_generator']['basis'], molecule=molecule)
    
    num_atoms = molecule.natom()

    with open(optimized_geometry_fn, 'w') as opt_geo_file:
        opt_geo_file.write(str(num_atoms) + "\n")
        opt_geo_file.write(energy_formatter.format(e) + "\n")
        
        for i in range(num_atoms):
            atom = molecule.symbol(i)
            
            if len(atom) > 1:
                atom = atom[:1] + atom[-1:].lower()
            

            x = filter_neg_zero(opt_formatter.format(molecule.x(i) / molecule.input_units_to_au()), 6) 
            y = filter_neg_zero(opt_formatter.format(molecule.y(i) / molecule.input_units_to_au()), 6)
            z = filter_neg_zero(opt_formatter.format(molecule.z(i) / molecule.input_units_to_au()), 6)
            
            opt_geo_file.write(atom + "\t" + x + " " + y + " " + z + "\n")
    
    print("")
    
    return molecule
    
def frequencies(molecule, config, normal_modes_fn):
    print("Determining normal modes and running frequency analysis...")

    norm_formatter = "{:> 12.4f}"
    freq_formatter = "{:.2f}"
    mass_formatter = "{:.6f}"
    
    total_energy, wavefunc = psi4.frequency(config['config_generator']['method'] + '/' + config['config_generator']['basis'], molecule=molecule, return_wfn=True)
    
    vib_info_raw = wavefunc.frequency_analysis
    vib_info = psi4.qcdb.vib.filter_nonvib(vib_info_raw)

    normal_out = ""
    
    frequencies = wavefunc.frequencies().to_array()
    reduced_masses = vib_info["mu"].data
    
    normal_modes = vib_info['x'].data

    num_modes = len(frequencies)
    num_atoms = molecule.natom()

    for i in range(1, 1 + num_modes):
        index = i - 1
        
        normal_out += "normal mode: " + str(i) + "\n"
        normal_out += freq_formatter.format(frequencies[index]) + "\n"
        normal_out += "red_mass = " + mass_formatter.format(reduced_masses[index]) + "\n"
        
        for atom in range(num_atoms):
            normal_out += filter_neg_zero(norm_formatter.format(float(normal_modes[3 * atom + 0][index])), 4) + "\t"
            normal_out += filter_neg_zero(norm_formatter.format(float(normal_modes[3 * atom + 1][index])), 4) + "\t"
            normal_out += filter_neg_zero(norm_formatter.format(float(normal_modes[3 * atom + 2][index])), 4) + "\n"
        
        normal_out += "\n"

    with open(normal_modes_fn, 'w') as norm_file:
        norm_file.write(normal_out)

    #print(vib_info)

    print("Normal mode/frequency analysis complete.")
    print("")

    return 3 * num_atoms - num_modes
