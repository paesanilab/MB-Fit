try:
    import psi4
    have_psi4 = True
except:
    have_psi4 = False

nthreads = 2

# MB-pol optimized water geometry
mol = """
 O      -6.738646395e-02  -1.491617397e+00  -1.095035971e-10
 H       8.141855337e-01  -1.866159045e+00  -2.044214632e-10
 H       7.436878209e-02  -5.443285903e-01   4.259078718e-11
"""

print(mol)

method = "BLYP"
basis = "cc-pvdz"
model = "{}/{}".format(method,basis)

if have_psi4:
    psi4.core.set_output_file("psi4.out", False)
    psi4.set_memory("1GB")
    psi4.geometry(mol)
    psi4.set_num_threads(nthreads)
    energy = psi4.energy(model)
    print("E({}) = {}".format(model, energy))
else:
    print("psi4 not available")
    
