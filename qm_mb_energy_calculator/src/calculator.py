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

try:
    from TensorMol import *
    os.environ["CUDA_VISIBLE_DEVICES"]="CPU" # for TensorMol
except ImportError:
    pass

import database

def sym_to_num(symbol):
    """
    Converts symbols into numbers for TensorMol input
    """
    if symbol == "O":
        return 8
    elif symbol == "H":
        return 1

def calc_energy(molecule, fragment_indicies, config):
    """
    Compute the energy using a model requested by the user
    """

# This section needs more thinking
# Now, instead of passing a string, it passes in an entire molecule 
    # Either create or connect to a database if a name is given
    connect, cursor = database(db)
    
    # Query once
    data = query(cursor, "Molecules", (molecule.to_xyz()))

    # Finish if we do not need to update
    if not config["database"]["update"] and data:
        database.finalize(connect)
        return some_property_of_data

    # If we have reached here, we intend to update data or it is not found
    model = config["driver"]["model"]
    
    if model == "psi4":
        psi4.core.set_output_file("/dev/null", False)
        psi4.set_memory(config["psi4"]["memory"])
        energy = calc_psi4_energy(molecule, fragment_indicies, config)
    
    elif model == "TensorMol":
        energy = TensorMol_convert_str(molecule, fragment_indicies, config)
  
    elif model == "qchem":
        energy = calc_qchem_energy(molecule, fragment_indicies, config)
    else:
        print("No such model exists!")
        return 0

    # If there was an existing entry, update
    if data:
        database.update(cursor, "Molecules", something)
    # Else, insert new entry
    else:
        database.insert(cursor, "Molecules", something)
    database.finalize(connect)
    return energy

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

def calc_psi4_energy(molecule, fragment_indicies, config):
    """
    Use PSI4 to compute the energy of a fragment
    """
    psi4_string = molecule.to_xyz(fragment_indicies)
    psi4_mol = psi4.core.Molecule.create_molecule_from_string(psi4_string)
    psi4_mol.update_geometry()

    psi4.set_num_threads(config["psi4"].getint("threads"))
    energy = psi4.energy("{}/{}".format(config["psi4"]["method"],
          config["psi4"]["basis"]), molecule=psi4_mol)
    
    #print(energy)
    return energy

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

    #prepare Q-chem input fil from frag string

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
