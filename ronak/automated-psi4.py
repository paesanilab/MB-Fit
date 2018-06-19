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


optimized_geometry_path = "inputs/optimized.xyz"
normal_modes_path = "inputs/normal_modes.dat"
output_path = "temp/psi4.out"

opt_formatter =  "{:< 16.10f}"
norm_formatter = "{:> 15.8f}"

if not has_psi4:
	print("ERROR: The psi4 module could not be imported.")
else:
	# Setting global psi4 variables
	memory = "1GB"
	num_threads = 2

	psi4.core.set_output_file(output_path, False)
	psi4.set_memory("1GB")
	psi4.set_num_threads(num_threads)

	# Defining the molecule
	mol = psi4.geometry("""
  		O 0.0 0.0 0.0
  		H 0.0 0.0 1.0
  		H 0.0 1.0 0.0
	""")

	# Defining the computations models for steps 1 and 2
	opt_method = "BLYP"
	opt_basis = "cc-pvdz"
	opt_model = opt_method + "/" + opt_basis
	
	freq_method = "BLYP"
	freq_basis = "cc-pvdz"
	freq_model = opt_method + "/" + opt_basis
	
	# Step 1
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
	total_energy, wavefunc = psi4.frequency(freq_model, molecule=mol, return_wfn=True)
	
#	In theory, this should get a vector containing python objects of the normal modes. Right now, it's just returning None
#
#	freqs = wavefunc.frequencies()
#	normal_modes = wavefunc.normalmodes()

	# Workaround that parses the output file
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

					normal_out += "normal mode: " + str(nor_mode_num) + "\n" + components[1] + "\nred_mass = -1.0\n"

				elif parse_lines > 0:
					parse_lines -= 1
					
					if "Force constant:" not in line and "mass" not in line:
						components = str.split(line)
						normal_out += norm_formatter.format(float(components[1])) + "\t" + norm_formatter.format(float(components[2])) + "\t" + norm_formatter.format(float(components[3])) + "\n"

				elif parse_lines == 0:
					normal_out += "\n"
					
	with open(normal_modes_path, 'w') as norm_file:
		norm_file.write(normal_out)


			
