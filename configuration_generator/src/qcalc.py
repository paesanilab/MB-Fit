# qcalc.py
#
# Calculator that uses accepts calls from nmcgen and calls upon the requested quantum chemistry code (e.g. psi4, qchem, etc) to carry out the specified calculation.

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../qm_mb_energy_calculator/src")
import io
import molecule_parser
from molecule import Molecule

try:
    import psi4
except ImportError:
    has_psi4 = False
else:
    has_psi4 = True

def init(config, log_name):
    verify_program(config)
    
    if config['config_generator']['code'] == 'psi4':
        psi4_helper.init(config, log_name)

def optimize(molecule, settings):
    
    if settings.get("config_generator", "code") == "psi4":
        return optimize_psi4(molecule, settings)
    
    elif config['config_generator']['code'] == 'qchem':
        energy, qchem_mol = qchem_helper.optimize(molecule, filenames, config)
        molecule = qchem_helper.read_qchem_mol(qchem_mol, molecule.num_atoms)
        return molecule, energy

def optimize_psi4(molecule, settings):

    psi4_string = "{}\n{} {}".format(molecule.to_xyz(), molecule.get_charge(), molecule.get_spin_multiplicity())

    try:
        psi4_mol = psi4.geometry(psi4_string)
    except RuntimeError as e:
        raise

    e = psi4.optimize(settings.get('config_generator', 'method') + '/' + settings.get('config_generator', 'basis'), molecule=psi4_mol)

    return Molecule().read_psi4_string(psi4_mol.save_string_xyz(), "molecule"), e

def frequencies(molecule, filenames, config):

    verify_program(config)

    if config['config_generator']['code'] == 'psi4':   
        psi4_molecule = psi4_helper.psi4_mol(molecule, config['molecule']['charges'], config['molecule']['spins'])
        return psi4_helper.frequencies(psi4_molecule, config)

    elif config['config_generator']['code'] == 'qchem':  
        return qchem_helper.frequencies(molecule, filenames, config)

