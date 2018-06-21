# automated-psi4.py
# 
# Automation of the following steps of the workflow:
#	1. geometry optimization (using psi4)
#	2. generation of normal modes/frequency calculation (using psi4)
#	3. configuration generation (using Sandra's code)
# Based on the user input of the initial geometry and the calculation model
#
# @author Ronak

try:
	import psi4
	has_psi4 = True
except:
	has_psi4 = False

import os
import sys

def treat_imaginary(contents):
	string = str(contents)
	
	result = string.replace("(","")
	result = result.replace(")","")
	
	plus_loc = result.index("+")
	result = result[:plus_loc]
	
	return result
		

optimized_geometry_path = "inputs/optimized.xyz"
normal_modes_path = "inputs/normal_modes.dat"

gcn_input_path = "gcn/conf_gen_input.txt"
gcn_output_path = "gcn/conf_gen_output.txt"

norm_config_path = "gcn/conf.xyz"

gcn_executor_path = "gcn/executor.sh"

# psi4 Global Parameters
memory = "1GB"
num_threads = 2
output_path = "temp/psi4.out"

opt_formatter =  "{:< 16.10f}"
norm_formatter = "{:> 15.8f}"

if len(sys.argv) == 2:
	config_passed = True
else:
	config_passed = False

if config_passed:
	# Parse the configuration file
	with open(sys.argv[1], 'r') as config_file:
		params = {"input_geo": 0, "linear_geo": 0, "charge": 0, "multiplicity": 0, "random": 0, "num_configs": 0, "geometric": 0, "linear": 0, "method": 0,"basis": 0}
		
		parse_fail = False
		
		try:
			for line in config_file:
				parts = line.split()
			
				if "input_geometry_path:" in line:
					with open(parts[1], 'r') as input_geo_file:
						input_geo = input_geo_file.read()
					params["input_geo"] += 1
					
				elif "linear_geo:" in line:
					linear_geo = parts[1] == "true"
					params["linear_geo"] += 1
					
				elif "charge:" in line:
					charge = parts[1]
					params["charge"] += 1
					
				elif "multiplicity:" in line:
					multiplicity = parts[1]
					params["multiplicity"] += 1
					
				elif "random:" in line:
					random = parts[1]
					params["random"] += 1
					
				elif "num_configs:" in line:
					num_configs = int(str(parts[1]))
					params["num_configs"] += 1
					
				elif "geometric:" in line:
					geometric = "." + parts[1] + "."
					params["geometric"] += 1
					
				elif "linear:" in line:
					linear = "." + parts[1] + "."
					params["linear"] += 1
					
				elif "method:" in line:
					method = parts[1]
					params["method"] += 1
					
				elif "basis:" in line:
					basis = parts[1]
					params["basis"] += 1
					
			for key in params:
				if params[key] != 1:
					parse_fail = True
		except:
			parse_fail = True

if not config_passed:
	print("Usage: automated-psi4.py [config.txt]")
elif not has_psi4:
	print("ERROR: The psi4 module could not be imported.")
elif parse_fail:
	print("ERROR: The configuration file does not conform to the input format.")
else:
	psi4.core.set_output_file(output_path, False)
	psi4.set_memory(memory)
	psi4.set_num_threads(num_threads)

	# Defining the molecule
	geo = input_geo + "\n" + charge + " " + multiplicity
	
	mol = psi4.geometry(geo)

	# Defining the computation models for steps 1 and 2
	opt_model = method + "/" + basis
	freq_model = method + "/" + basis
	
	if linear_geo:
		dim_null = 5
	else:
		dim_null = 6
	
	# Step 1
	print("Optimizing geometry...")
	e = psi4.optimize(opt_model, molecule=mol)
	
	num_atoms = mol.natom()

	with open(optimized_geometry_path, 'w') as opt_geo_file:
		opt_geo_file.write(str(num_atoms) + "\n")
		opt_geo_file.write(str(e) + "\n")
		
		for i in range(num_atoms):
			atom = mol.symbol(i)

			x = opt_formatter.format(mol.x(i))
			y = opt_formatter.format(mol.y(i))
			z = opt_formatter.format(mol.z(i))
			
			opt_geo_file.write(atom + "\t" + x + " " + y + " " + z + "\n")
	
	# Step 2
	print("")
	print("Determining normal modes and running frequency analysis...")
	total_energy, wavefunc = psi4.frequency(freq_model, molecule=mol, return_wfn=True)
	
	# Using the nightly build, we can actually get the normal modes in python
	vib_info = wavefunc.frequency_analysis
	
	num_modes = 3*num_atoms - dim_null
	
	additive = len(vib_info["omega"].data) - num_modes - 1

	normal_out = ""

	for i in range(1, 1 + num_modes):
		index = additive + i
		
		normal_out += "normal mode: " + str(i) + "\n"
		normal_out += treat_imaginary(vib_info["omega"].data[index]) + "\n"
		normal_out += "red_mass = " + str(vib_info["mu"].data[index]) + "\n"
		
		coords = vib_info["q"].data[index]
		for i in range(num_atoms):
			normal_out += norm_formatter.format(float(coords[3 * i])) + "\t" + norm_formatter.format(float(coords[3 * i + 1])) + "\t" + norm_formatter.format(float(coords[3 * i + 2])) + "\n"
		
		normal_out += "\n"

	ignore = """
	# Otherwise, we're going to need a stable v1.1 workaround that parses the output file
	normal_out = ""

	with open(output_path, 'r') as out_file:
		parseing = False
		nor_mode_num = 0
		parse_lines = -1

		for line in out_file:
			if "Normal Modes (non-mass-weighted)." in line:
				parseing = True

			if "-------------------------------------------------------------" in line:
				parseing = False
			
			if parseing:
				if "Frequency:" in line:
					components = str.split(line)

					nor_mode_num += 1
					parse_lines = num_atoms + 2

					normal_out += "normal mode:  " + str(nor_mode_num) + "\n" + components[1] + "\nred_mass = -1.0\n"

				elif parse_lines > 0:
					parse_lines -= 1
					
					if "Force constant:" not in line and "mass" not in line:
						components = str.split(line)
						normal_out += norm_formatter.format(float(components[1])) + "\t" + norm_formatter.format(float(components[2])) + "\t" + norm_formatter.format(float(components[3])) + "\n"

				elif parse_lines == 0:
					normal_out += "\n"
"""

	with open(normal_modes_path, 'w') as norm_file:
		norm_file.write(normal_out)
	print("Normal mode/frequency analysis complete.")
	print("")

	# Step 3
	print("Running normal distribution configuration generator...")
	with open(gcn_input_path, 'w') as input_file:
		input_file.write("\'" + optimized_geometry_path + "\'\n") 			# optimized geometry
		input_file.write("\'" + normal_modes_path + "\'\n") 				# normal modes
		input_file.write(str(3 * num_atoms) + " " + str(dim_null) + "\n")	# dim; dimnull
		input_file.write(random + " " + str(num_configs) + "\n")			# random method; nconfigs
		input_file.write("\'" + norm_config_path + "\'\n")					# output one-body configurations
		input_file.write(geometric + " " + linear + "\n")					# geometric, linear
		input_file.write(".true.")											# verbose

	with open(gcn_executor_path, 'w') as executor_file:
		executor_file.write("#!/bin/bash\n\n")
		executor_file.write("gcn/src/generate_configs_normdistrbn < " + gcn_input_path + " > " + gcn_output_path + "\n")
		executor_file.write("exit 0")
		
	os.system("bash " + gcn_executor_path)
	print("Normal Distribution Configuration generation complete.")

			
