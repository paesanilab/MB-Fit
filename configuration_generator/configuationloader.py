# configuationloader
# 
# Parses the configation text file, using the configparser module
#
# @author Ronak

from configparser import SafeConfigParser

def parse_config(path):
	parser = SafeConfigParser()
	parser.read(path)

	parse_fail = False
	
	name = None
	input_geo = ""
	linear_geo = None
	charge = None
	multiplicity = None
	random = None
	num_configs = None
	geometric = None
	linear = None
	method = None
	basis = None
	
	
	try:
	
		input_geo_path = parser.get("file_input", "input_geometry_path")
				
		with open(input_geo_path, 'r') as input_file:
			for i, line in enumerate(input_file):
				if i > 1:
					input_geo += line	
		
		if input_geo_path.endswith('.xyz'):
			name = input_geo_path[:-4]
		else:
			parse_fail = True
		
		linear_geo = parser.get("molecule", "linear_geo") == "true"
		charge = parser.get("molecule", "charge")
		multiplicity = parser.get("molecule", "multiplicity")
		
		random = parser.get("config_generator", "random")
		num_configs = parser.get("config_generator", "num_configs")
		geometric = parser.get("config_generator", "geometric")
		linear = parser.get("config_generator", "linear")
		
		method = parser.get("psi4", "method")
		basis = parser.get("psi4", "basis")
	except:
		parse_fail = True
			
	return parse_fail, name, input_geo, linear_geo, charge, multiplicity, random, num_configs, geometric, linear, method, basis 
