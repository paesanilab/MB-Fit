"""
Test script for minimal psi4 calculation
"""
import sys
sys.path.insert(0, "../src")

from molecule import Fragment, Molecule
import psi4

minimal_frag = Fragment()

minimal_frag.natoms = 3
minimal_frag.symbols = ['O','H','H']
minimal_frag.coordinates =[
[-0.0005900753,0.0000000000,-0.0004526740],
[0.9256102457,0.0000000000,-0.2481132352],
[-0.0000930743,0.0000000000,0.9582869092]
]

minimal_mol = Molecule()
minimal_mol.natoms = 3
minimal_mol.fragments = [minimal_frag]

psi4.core.set_output_file("/dev/null", False)
psi4.set_memory("1GB")

psi4_mol =  psi4.core.Molecule.create_molecule_from_string(
    minimal_mol.mol_comb_str([0]))
psi4_mol.update_geometry()
psi4.set_num_threads(6)
energy = psi4.energy("HF/STO-3G", molecule=psi4_mol)
    
print(energy)

