# qcalc.py
#
# Calculator that uses accepts calls from nmcgen and calls upon the requested quantum chemistry code (e.g. psi4, qchem, etc) to carry out the specified calculation.

import psi4_helper
import psi4

from molecule import Molecule

supported_programs = ['psi4'] #, 'qchem']

def init(config, log_name):
    if config['program']['code'] not in supported_programs:
        raise ValueError(program + " is not the name of a supported program")
        
    if config['program']['code'] == "psi4":
        psi4.core.set_output_file(log_name + ".log", False)
        psi4.set_memory(config['program']['memory'])
        psi4.set_num_threads(int(config['program']['num_threads']))

def verify_program(config):
    if config['program']['code'] not in supported_programs:
        raise ValueError(program + " is not the name of a supported program")
        
def optimize(molecule, config):
    verify_program(config)
    
    if config['program']['code'] == "psi4":
        psi4_molecule, energy = psi4_helper.optimize(molecule.psi4(config), config)
         
        molecule.read_psi4_mol(psi4_molecule)
        
        return molecule, energy
        
def frequencies(molecule, config):
    verify_program(config)
    
    if config['program']['code'] == "psi4":
        return psi4_helper.frequencies(molecule.psi4(config), config)
