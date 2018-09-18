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

import os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import settings_reader
import geometry_optimizer
import normal_modes_generator
import configuration_generator

def nmcgen(settings_path, unopt_geo_path, opt_geo_path, normal_modes_path, configs_path):
    geometry_optimizer.optimize_geometry(settings_path, unopt_geo_path, opt_geo_path)
    dim_null = normal_modes_generator.generate_normal_modes(settings_path, opt_geo_path, normal_modes_path)
    configuration_generator.generate_configurations(settings_path, opt_geo, normal_modes, dim_null, configs_path)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python nmcgen.py <settings path> <unopt geo path> <opt geo path> <normal modes path> <configs path>")
        exit(1)
    nmcgen(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
