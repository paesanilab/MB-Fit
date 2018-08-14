"""
A module for the different models for calculating the 
energy of a set of atoms (a fragment)
"""
import os, sys, subprocess
import numpy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")

import settings_reader
import constants

from exceptions import LibraryNotAvailableError, LibraryCallError, NoSuchLibraryError, ConfigMissingSectionError, ConfigMissingPropertyError

from psi4.driver.qcdb.exceptions import QcdbException


has_psi4 = True

try:
    import psi4
except ImportError:
    has_psi4 = False

# this is commented because we are not making TensorMol a priority
'''
try:
    from TensorMol import *
    os.environ["CUDA_VISIBLE_DEVICES"]="CPU" # for TensorMol
except ImportError:
    pass
'''

def calculate_energy(molecule, fragment_indicies, model, cp, settings):
    """
    Calculates the energy of a subset of the fragments of a molecule.

    Args:
        molecule    - the Molecule object to calculate the energy of
        fragment_indicies - list of indicies of fragments to include in the calculation
        model       - the model to use for this calculation, should be specified as method/basis
        cp          - whether to use counterpoise correction, True means that fragments not included in
                    fragment_indicies will be added as ghost atoms to the calculation
        settings    - .ini file with settings information

    Returns:
        The calculated energy
    """

        # check which library to use to calculate the energy
    code = settings.get("energy_calculator", "code")
    
    if code == "psi4":
        return calc_psi4_energy(molecule, fragment_indicies, model, cp, settings)
    
    elif code == "TensorMol":
        return TensorMol_convert_str(molecule, fragment_indicies, settings)
  
    elif code == "qchem":
        return calc_qchem_energy(molecule, fragment_indicies, model, cp, settings)
    
    # if the user has specified a code that does not exist, raise an Exception
    else:
        raise NoSuchLibraryError(code)


def calc_psi4_energy(molecule, fragment_indicies, model, cp, settings):
    """
    Calculates the energy of a subset of the fragments of a molecule with psi4

    Args:
        molecule    - the Molecule object to calculate the energy of
        fragment_indicies - list of indicies of fragments to include in the calculation
        model       - the model to use for this calculation, should be specified as method/basis
        cp          - whether to use counterpoise correction, True means that fragments not included in
                    fragment_indicies will be added as ghost atoms to the calculation
        settings    - .ini file with settings information

    Returns:
        The calculated energy
    """
    
    # if psi4 is not loaded, raise an exception
    if not has_psi4:
        raise LibraryNotAvailableError("psi4")

    # file to write logs from psi4 calculation
    log_file = settings.get("files", "log_path") + "/calculations/" + model + "/" + str(cp) + "/" + molecule.get_SHA1()[:8] + "frags:" + fragments_to_energy_key(fragment_indicies) + ".out"
    
    # create log file's directory if it does not already exist
    if not os.path.exists(os.path.dirname(log_file)):
        os.makedirs(os.path.dirname(log_file))
    
    # set the log file
    psi4.core.set_output_file(log_file, False)

    # Creats the psi4 input string of the molecule by combining the xyz file output with an additional line containing charge and spin multiplicity
    psi4_string = molecule.to_xyz(fragment_indicies, cp) + "\n" + str(molecule.get_charge(fragment_indicies)) + " " + str(molecule.get_spin_multiplicity(fragment_indicies))

    try:
        # Constructs the psi4 Molecule from the string representation by making a call to a psi4 library function
        psi4_mol = psi4.geometry(psi4_string)
    except RuntimeError as e:
        raise LibraryCallError("psi4", "geometry", str(e))

    # Set the number of threads to use to compute based on settings
    psi4.set_num_threads(settings.getint("psi4", "num_threads"))

    # Set the amount of memory to use based on settings
    psi4.set_memory(settings.get("psi4", "memory"))
    
    try:
        # Perform library call to calculate energy of molecule
        return psi4.energy(model, molecule=psi4_mol)
    except QcdbException as e:
        raise LibraryCallError("psi4", "energy", str(e))

def TensorMol_convert_str(molecule, fragment_indicies, settings):
    """
    Prepares input data 
    """
    tensorMol_string = molecule.to_xyz(fragment_indicies)
    xyz_list = tensorMol_string.split()
    atoms = numpy.array([])
    coords = numpy.array([])
    for i in range(int(len(xyz_list)/4)):
        atoms = numpy.append(atoms, constants.symbol_to_number(xyz_list[i*4]))
        coords = numpy.append(coords, xyz_list[i*4+1:(i+1)*4])
    coords = coords.reshape(atoms.size, 3)
    return calc_TensorMol_energy(atoms, coords, settings)

def calc_TensorMol_energy(atoms, coords, settings):
    """
    Sets up calculation for TensorMol
    """
    molecule = Mol(atoms, coords)
    molset = MSet()
    molset.mols.append(molecule)
    manager = GetWaterNetwork(molset)
    return En(molecule, coords, manager)

def calc_qchem_energy(molecule, fragment_indicies, model, cp, settings):
    """
    Calculates the energy of a subset of the fragments of a molecule with q-chem

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

    # Check if the user has qchem isntalled
    if subprocess.call("which qchem", shell=True) != 0:
        raise LibraryNotAvailableError("qchem")

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
    qchem_input += molecule.to_xyz(fragment_indicies, cp) + "\n"

    qchem_input += "$end\n"

    # Q-chem settings
    qchem_input += "$rem\n"

    qchem_input += "jobtype " + "sp" + "\n"
    qchem_input += "method " + model.split('/')[0] + "\n"
    qchem_input += "basis " + model.split('/')[1] + "\n"

    qchem_input += "ecp " + settings.get("energy_calculator", "ecp") + "\n"

    qchem_input += "$end"

    # file to write qchem input in
    log_file_in = settings.get("files", "log_path") + "/calculations/" + model + "/" + str(cp) + "/" + molecule.get_SHA1()[:8] + "frags:" + fragments_to_energy_key(fragment_indicies) + ".in"
    
    # make sure directory exists prior to writing to file
    if not os.path.exists(os.path.dirname(log_file_in)):
        os.makedirs(os.path.dirname(log_file_in))
    f = open(log_file_in, "w")
    f.write(qchem_input)
    f.close()
    
    # file to write qchem output in
    log_file_out = settings.get("files", "log_path") + "/calculations/" + model + "/" + str(cp) + "/" + molecule.get_SHA1()[:8] + "frags:" + fragments_to_energy_key(fragment_indicies) + ".out"

    # get number of threads
    num_threads == settings.get("qchem", "num_threads")
    

    # perform system call to run qchem
    # want to save log stuff rather than put in /dev/null?
    if subprocess.call("qchem " +  log_file_in + " " + log_file_out + " -nt " + num_threads + " >> /dev/null", shell=True) != 0:
        raise LibraryCallError("qchem", "energy calculation", "process returned non-zero exit code")
    
    # find the energy inside the qchem file output
    # the with open as syntax automatically closes the file
    with open(log_file_out) as qchem_out:
        for line in qchem_out:
            if line.find("Total energy in the final basis set = ") != -1:
                return float(line[line.find("Total energy in the final basis set = ") + 39:])

    # if no line with "Total energy in the final basis set = " is found, raise an exception
    raise LibraryCallError("qchem", "energy calculation", "process returned file of incorrect format")

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

def fragments_to_energy_key(fragments):
    """
    Generate a simple string from a list of fragments. [1,2,3] -> E123, [0] -> E0
    For use in naming of log files

    Args:
        fragments   - list of indicies of fragments

    Returns:
        string of fromat E123 where 123 are the numbers specified by the fragments list
    """

    # initialize the key
    key = "E"

    # loop thru every fragment index in the list
    for fragment in fragments:
        key += str(fragment)

    return key
