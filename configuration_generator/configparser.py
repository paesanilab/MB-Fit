# configparser
# 
# Parses the configation text file, using a format that I made up. Will be replaced by the actual widespread library configparser
#
# @author Ronak

# Parses the config.txt file, uses a methodology that I invented.
# TODO: Implement configparser
def parse_config(path):
	parse_fail = False
	
	input_geo = None
	linear_geo = None
	charge = None
	multiplicity = None
	random = None
	num_configs = None
	geometric = None
	linear = None
	method = None
	basis = None
	
	with open(path, 'r') as config_file:
		params = {"input_geo": 0, "linear_geo": 0, "charge": 0, "multiplicity": 0, "random": 0, "num_configs": 0, "geometric": 0, "linear": 0, "method": 0, "basis": 0}
		
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
			
	return parse_fail, input_geo, linear_geo, charge, multiplicity, random, num_configs, geometric, linear, method, basis 
