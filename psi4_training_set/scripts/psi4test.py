import os, sys
import psi4

psi4.core.set_output_file('output.dat', False)
psi4.set_memory('500 MB')
h2o = psi4.geometry("""
O
H 1 0.96  
H 1 0.96 2 104.5  
""")

psi4.energy('scf/cc-pvdz')
