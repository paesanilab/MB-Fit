# psi4_helper
# 
# Allows for function calls to the Psi4 library in order to optimize the molecular geometry and to generate the normal modes/caluclate their frequencies.
#
# @author Ronak

import os
import sys

from molecule import Molecule
import numpy as np

try:
    import psi4
except:
    print("psi4 not found")


def init(config, log_name):
    psi4.core.set_output_file(log_name + ".log", False)
    psi4.set_memory(config['program']['memory'])
    psi4.set_num_threads(int(config['program']['num_threads']))

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
    
def read_psi4_mol(psi4_mol, au_conversion=1.0):
    return Molecule(psi4_mol.create_psi4_string_from_molecule(), au_conversion)

def psi4_mol(molecule, charge, multiplicity):
    geo = str(molecule) + "\n" + charge + " " + multiplicity
        
    return psi4.geometry(geo)
