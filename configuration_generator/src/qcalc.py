# qcalc.py
#
# Calculator that uses accepts calls from nmcgen and calls upon the requested quantum chemistry code (e.g. psi4, qchem, etc) to carry out the specified calculation.

import psi4_helper
from molecule import Molecule

supported_programs = ['psi4'] #, 'qchem']

def verify_program(config):
    if config['program']['code'] not in supported_programs:
        raise ValueError(program + " is not the name of a supported program")
        
def optimize(molecule, config, optimized_geometry_fn):
    verify_program(config)
    
    if config['program']['code'] == "psi4":
        psi4_molecule = psi4_helper.optimize(molecule.psi4(), config, optimized_geometry_fn)
        
        molecule.read_psi4_mol(psi4_molecule)
        
def frequencies(molecule, config, normal_modes_fn):
    verify_program(config)
    
    if config['program']['code'] == "psi4":
        return psi4_helper.frequencies(molecule.psi4(), config, normal_modes_fn)
