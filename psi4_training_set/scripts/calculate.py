import psi4

def calculate(mol, args):
    """ psi4 calculations file """

    # Some constant args from the config
    memory = args.memory
    method = args.method
    basis = args.basis
    charge = args.charge
    spin = args.spin
    convergence = args.convergence
    cycles = args.cycles
    threshold = args.threshold

    try:
        psi4.set_memory(memory)
    except:
        print("Invalid memory settings.")
        return None

    # Something calculation-related
    psi4_mol = psi4.core.Molecule.create_molecule_from_string(str(mol)) 

    psi4_mol.update_geometry()
   
    # Given the test input, we are calculating the same calculations 10 times.
    # There should be a way to only calculate this once.
    #storing the single-point electronic energy in a variable
    try:
        ref_energy = psi4.energy("{}/{}".format(method, basis), 
            molecule=psi4_mol)
    except:
        print("Method does not exist.")
        return None

    # Store the one-body energy into molecule
    mol.set_energy(ref_energy)
    print(mol.energy)
    
    # Return the energy for file to write out
    return ref_energy
