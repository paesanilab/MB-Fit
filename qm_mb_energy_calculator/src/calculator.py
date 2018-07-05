"""
A module for the different models for calculating the 
energy of a set of atoms (a fragment)
"""
try:
    import numpy
except:
    pass

from subprocess import call

try:
    import psi4
except ImportError:
    print("Error: failed to import psi4")
    pass

'''
try:
    from TensorMol import *
    os.environ["CUDA_VISIBLE_DEVICES"]="CPU" # for TensorMol
except ImportError:
    pass
'''

def sym_to_num(symbol):
    """
    Converts symbols into numbers for TensorMol input
    """
    if symbol == "O":
        return 8
    elif symbol == "H":
        return 1


def calculate_energy(molecule, fragment_indicies, model, cp, config):
    """
    Computes energy of the given fragments of a molecule.

    Uses only the fragments indicated by the fragment indices array.

    Parameters
    ----------
    arg1: Molecule
        The Molecule to calculate the energy of.
    arg2: int list
        Which fragments to take into account when calculating energy.
    arg3: str
        What model to use to calcuate the enrgy. Saved as method/basis
    arg4: boolean
        Whether to use Superposition Corrections for the basis set
    arg5: ConfigParser
        Contains settings, etc, that dictate how to calculate energy.

    Returns
    -------
    int
        Calculated Energy.
    """

    # check which library to use to calculate the energy
    code = config["energy_calculator"]["code"]
    
    if code == "psi4":
        # throw away output
        psi4.core.set_output_file("/dev/null", False)
        return calc_psi4_energy(molecule, fragment_indicies, model, cp, config)
    
    elif code == "TensorMol":
        return TensorMol_convert_str(molecule, fragment_indicies, config)
  
    elif code == "qchem":
        return calc_qchem_energy(molecule, fragment_indicies, config)
    else:
        raise ValueError("{} is not a valid code library".format(code))


def calc_psi4_energy(molecule, fragment_indicies, model, cp, config):
    """
    Uses psi4 library to compute energy of a molecule

    Uses only the fragments indicated by the fragment indices array.

    Parameters
    ----------
    arg1: Molecule
        The Molecule to calculate the energy of.
    arg2: int list
        Which fragments to take into account when calculating energy.
    arg3: str
        What model to use to calcuate the enrgy. Saved as method/basis
    arg4: boolean
        Whether to use Superposition Corrections for the basis set
    arg5: ConfigParser
        Contains settings, etc, that dictate how to calculate energy.

    Returns
    -------
    int
        Calculated Energy.
    """

    # Creats the psi4 input string of the molecule by combining the xyz file output with an additional line containing charge and spin multiplicity
    psi4_string = molecule.to_xyz(fragment_indicies) + "\n" + str(molecule.get_charge()) + " " + str(molecule.get_spin())

    # Constructs the psi4 Molecule from the string representation by making a call to a psi4 library function
    psi4_mol = psi4.core.Molecule.create_molecule_from_string(psi4_string)

    # Updates geometry of psi4 Molecule. !!! Unsure what this does !!!
    psi4_mol.update_geometry()

    # Set the number of threads to use to compute based on configuration
    psi4.set_num_threads(config["psi4"].getint("num_threads"))

    # Set the amount of memory to use based on configuration
    psi4.set_memory(config["psi4"]["memory"])

    # set the name of the bsse_type to calculate energy
    if cp == "True":
        cp_type = "cp"
    else:
        cp_type = "nocp"

    # Perform library call to calculate energy of molecule
    return psi4.energy(model, molecule=psi4_mol)

def TensorMol_convert_str(molecule, fragment_indicies, config):
    """
    Prepares input data 
    """
    tensorMol_string = molecule.to_xyz(fragment_indicies)
    xyz_list = tensorMol_string.split()
    atoms = numpy.array([])
    coords = numpy.array([])
    for i in range(int(len(xyz_list)/4)):
        atoms = numpy.append(atoms, sym_to_num(xyz_list[i*4]))
        coords = numpy.append(coords, xyz_list[i*4+1:(i+1)*4])
    coords = coords.reshape(atoms.size, 3)
    return calc_TensorMol_energy(atoms, coords, config)

def calc_TensorMol_energy(atoms, coords, config):
    """
    Sets up calculation for TensorMol
    """
    molecule = Mol(atoms, coords)
    molset = MSet()
    molset.mols.append(molecule)
    manager = GetWaterNetwork(molset)
    return En(molecule, coords, manager)

"""
Calculates energy using qchem
"""
def calc_qchem_energy(molecule, fragment_indicies, config):
    """
    Uses q-chem library to compute energy of a molecule

    Uses only the fragments indicated by the fragment indices array.

    Parameters
    ----------
    arg1: Molecule
        The Molecule to calculate the energy of.
    arg2: int list
        Which fragments to take into account when calculating energy.
    arg3: ConfigParser
        Contains settings, etc, that dictate how to calculate energy.

    Returns
    -------
    int
        Calculated Energy.
    """

    #prepare Q-chem input file from frag string

    input_file = "qchem_input.txt"
    output_file = "qchem_out.txt"
    
    # initialize qchem input string
    qchem_input = "";
    
    # molecule format
    qchem_input += "$molecule\n"

    # charge and spin multiplicity
    qchem_input += "0 1\n"

    # atoms in the molecule
    # might need to add whitespace before each line?
    qchem_input += molecule.to_xyz(fragment_indicies)

    qchem_input += "$end\n"

    # Q-chem settings
    qchem_input += "$rem\n"

    qchem_input += "jobtype " + "sp" + "\n"
    qchem_input += "method " + config["Qchem"]["method"] + "\n"
    qchem_input += "basis " + config["Qchem"]["basis"] + "\n"

    qchem_input += "$end"
 
    call(["mkdir", "qchem"]);
    
    f = open("qchem/qchem_in",'w');
    f.write(qchem_input)
    f.close();    
    
    # call(["qchem", "qchem/qchem_in", "qchem/qchem_out", "> /dev/null"]);
    return 0

# Water network data is required to be in the same directory under ./networks !!
def GetWaterNetwork(a):
    """
    TensorMol function that loads a pretrained water network
    """
    TreatedAtoms = a.AtomTypes()
    PARAMS["tf_prec"] = "tf.float64"
    PARAMS["NeuronType"] = "sigmoid_with_param"
    PARAMS["sigmoid_alpha"] = 100.0
    PARAMS["HiddenLayers"] = [500, 500, 500]
    PARAMS["EECutoff"] = 15.0
    PARAMS["EECutoffOn"] = 0
    PARAMS["Elu_Width"] = 4.6  # when elu is used EECutoffOn should always equal to 0
    PARAMS["EECutoffOff"] = 15.0
    PARAMS["DSFAlpha"] = 0.18
    PARAMS["AddEcc"] = True
    PARAMS["KeepProb"] = [1.0, 1.0, 1.0, 1.0]
    d = MolDigester(TreatedAtoms, name_="ANI1_Sym_Direct", OType_="EnergyAndDipole")
    tset = TensorMolData_BP_Direct_EE_WithEle(a, d, order_=1, num_indis_=1, type_="mol",  WithGrad_ = True)
    manager=TFMolManage("water_network",tset,False,"fc_sqdiff_BP_Direct_EE_ChargeEncode_Update_vdw_DSF_elu_Normalize_Dropout",False,False)
    return manager

def En(m, x_, manager):
    """
    Use TensorMol to compute the energy of a fragment
    """
    mtmp = Mol(m.atoms,x_)
    Etotal, Ebp, Ebp_atom, Ecc, Evdw, mol_dipole, atom_charge, gradient = manager.EvalBPDirectEEUpdateSingle(mtmp, PARAMS["AN1_r_Rc"], PARAMS["AN1_a_Rc"],PARAMS["EECutoffOff"], True)
    energy = Etotal[0]
    return energy

