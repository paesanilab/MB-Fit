# psi4_helper
# 
# Allows for function calls to the Psi4 library in order to optimize the molecular geometry and to generate the normal modes/caluclate their frequencies.
#
# @author Ronak

import psi4
import os
import sys

import numpy as np

def optimize(molecule, config):
    print("Optimizing geometry...")

    e = psi4.optimize(config['model']['method'] + '/' + config['model']['basis'], molecule=molecule)
    
    print("")
    
    return molecule, e
    
def frequencies(molecule, config):
    print("Determining normal modes and running frequency analysis...")
    
    total_energy, wavefunc = psi4.frequency(config['model']['method'] + '/' + config['model']['basis'], molecule=molecule, return_wfn=True)
    
    vib_info_raw = wavefunc.frequency_analysis
    vib_info = psi4.qcdb.vib.filter_nonvib(vib_info_raw)
    
    frequencies = wavefunc.frequencies().to_array()
    red_masses = vib_info["mu"].data
    
    normal_modes_raw = vib_info['x'].data
    normal_modes = []
    
    num_modes = len(frequencies)
    num_atoms = molecule.natom()

    for i in range(num_modes):
        geometry = np.zeros((num_atoms, 3))
        
        for atom in range(num_atoms):
            geometry[atom][0] = float(normal_modes_raw[3 * atom + 0][i])
            geometry[atom][1] = float(normal_modes_raw[3 * atom + 1][i])
            geometry[atom][2] = float(normal_modes_raw[3 * atom + 2][i])
        
        normal_modes.append(geometry)

    print("Normal mode/frequency analysis complete.")
    print("")
    
    return normal_modes, frequencies, red_masses
