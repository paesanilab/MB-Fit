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

# Main Functions
def optimize_geometry(molecule, model, optimized_geometry_path):
	print("Optimizing geometry...")

	opt_formatter =  "{:< 16.10f}"

	e = psi4.optimize(model, molecule=molecule)
	
	num_atoms = molecule.natom()

	with open(optimized_geometry_path, 'w') as opt_geo_file:
		opt_geo_file.write(str(num_atoms) + "\n")
		opt_geo_file.write(str(e) + "\n")
		
		for i in range(num_atoms):
			atom = molecule.symbol(i)

			x = opt_formatter.format(molecule.x(i))
			y = opt_formatter.format(molecule.y(i))
			z = opt_formatter.format(molecule.z(i))
			
			opt_geo_file.write(atom + "\t" + x + " " + y + " " + z + "\n")
			
	print("")
	
def frequency_calculations(molecule, model, dim_null, normal_modes_path):
	print("Determining normal modes and running frequency analysis...")

	norm_formatter = "{:> 15.8f}"

	total_energy, wavefunc = psi4.frequency(model, molecule=molecule, return_wfn=True)
	
	vib_info_raw = wavefunc.frequency_analysis
	vib_info = psi4.qcdb.vib.filter_nonvib(vib_info_raw)
	
	num_atoms = molecule.natom()
	num_modes = 3*num_atoms - dim_null

	normal_out = ""
	
	frequencies = vib_info["omega"].data
	reduced_masses = vib_info["mu"].data
	
	normal_modes = vib_info['x'].data

	for i in range(1, 1 + num_modes):
		index = i - 1
		
		normal_out += "normal mode: " + str(i) + "\n"
		normal_out += treat_imaginary(frequencies[index]) + "\n"
		normal_out += "red_mass = " + str(reduced_masses[index]) + "\n"
		
		for atom in range(num_atoms):
			normal_out += norm_formatter.format(float(normal_modes[3 * atom + 0][index])) + "\t"
			normal_out += norm_formatter.format(float(normal_modes[3 * atom + 1][index])) + "\t"
			normal_out += norm_formatter.format(float(normal_modes[3 * atom + 2][index])) + "\n"
		
		normal_out += "\n"


	with open(normal_modes_path, 'w') as norm_file:
		norm_file.write(normal_out)

	#print(vib_info)

	print("Normal mode/frequency analysis complete.")
	print("")
