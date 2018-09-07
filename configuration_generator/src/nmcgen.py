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

import os
import sys
import shutil

from configparser import ConfigParser

from molecule import Molecule
import qcalc

import gcn_runner
import config_loader

import output_writer

config = config_loader.load()
filenames, log_name = config_loader.process_files(config)
  
qcalc.init(config, log_name)

def nmcgen(settings_path, unopt_geo, opt_geo, normal_modes, configs):
    settings = settings_reader.SettingsReader(settings_path)
    geometry_optimizer.optimize_geometry(settings_path, unopt_geo, opt_geo)
    dim_null = normal_modes_generator.generate_normal_modes(settings_path, opt_geo, normal_modes)
    configuration_generator.generate_configs(settings_path, opt_geo, normal_modes, dim_null, configs)

    normal_modes_generator.generate_normalModes(settings_path, geometry)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python nmcgen.py <settings path> <unopt geo> <opt geo> <normal modes> <configs>")
        exit(1)
    nmcgen(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
