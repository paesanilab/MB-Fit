# nmcgen.py
# 
# Normal Mode Configuration Generator
#
# Automation of the following steps of the workflow:
#	1. geometry optimization
#	2. generation of normal modes/frequency calculation
#	3. configuration generation
#
# The calculations for steps one and two are completed using a quantum chemistry code, psi4 by default. Step three is completed using Sandra's GCN code.
#
# @author Ronak

import sys

from potential_fitting import configurations

def nmcgen(settings_path, unopt_geo_path, opt_geo_path, normal_modes_path, configs_path):
    configurations.optimize_geometry(settings_path, unopt_geo_path, opt_geo_path)
    dim_null = configurations.generate_normal_modes(settings_path, opt_geo_path, normal_modes_path)
    configurations.generate_1b_configurations(settings_path, opt_geo_path, normal_modes_path, dim_null, configs_path)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python nmcgen.py <settings path> <unopt geo path> <opt geo path> <normal modes path> <configs path>")
        exit(1)
    nmcgen(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
