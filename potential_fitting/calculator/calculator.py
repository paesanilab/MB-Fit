# external package imports
import os, subprocess

# absolute module imports
from potential_fitting.utils import SettingsReader, constants, files, system
from potential_fitting.exceptions import LibraryNotAvailableError, LibraryCallError, NoSuchLibraryError, \
        ConfigMissingSectionError, ConfigMissingPropertyError, CommandExecutionError

has_psi4 = True

try:
    from psi4.driver.qcdb.exceptions import QcdbException
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
        molecule            - The Molecule object to calculate the energy of.
        fragment_indicies   - List of indicies of fragments to include in the calculation.
        model               - The model to use for this calculation, should be specified as method/basis.
        cp                  - Whether to use counterpoise correction, True means that fragments not included in
                fragment_indicies will be added as ghost atoms to the calculation.
        settings            - SettingsReader object with all relevent settings.

    Returns:
        The calculated energy of the subset of fragments in this molecule.
    """

        # check which library to use to calculate the energy
    code = settings.get("energy_calculator", "code")
    
    if code == "psi4":
        return calc_psi4_energy(molecule, fragment_indicies, model, cp, settings)
    '''
    elif code == "TensorMol":
        return TensorMol_convert_str(molecule, fragment_indicies, settings)
    '''
    elif code == "qchem":
        return calc_qchem_energy(molecule, fragment_indicies, model, cp, settings)
    
    # if the user has specified a code that does not exist, raise an Exception
    else:
        raise NoSuchLibraryError(code)


def calc_psi4_energy(molecule, fragment_indicies, model, cp, settings):
    """
    Calculates the energy of a subset of the fragments of a molecule with psi4.

    Args:
        molecule            - The Molecule object to calculate the energy of.
        fragment_indicies   - List of indicies of fragments to include in the calculation.
        model               - The model to use for this calculation, should be specified as method/basis.
        cp                  - Whether to use counterpoise correction, True means that fragments not included in
                    fragment_indicies will be added as ghost atoms to the calculation.
        settings            - SettingsReader object with all relevent settings.

    Returns:
        The calculated energy of the subset of fragments in this molecule.
    """

    method, basis = model.split("/")
    
    # if psi4 is not loaded, raise an exception
    if not has_psi4:
        raise LibraryNotAvailableError("psi4")

    # file to write logs from psi4 calculation
    log_path = files.get_energy_log_path(settings.get("files", "log_path"), molecule, method, basis, cp, "out")
    
    # set the log file
    psi4.core.set_output_file(log_path, False)

    # Creats the psi4 input string of the molecule by combining the xyz file output with an additional line containing charge and spin multiplicity
    psi4_string = "{}\n{} {}".format(molecule.to_xyz(fragment_indicies, cp), molecule.get_charge(fragment_indicies),
            molecule.get_spin_multiplicity(fragment_indicies))

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

'''
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
'''

def calc_qchem_energy(molecule, fragment_indicies, model, cp, settings):
    """
    Calculates the energy of a subset of the fragments of a molecule with q-chem.

    Args:
        molecule            - The Molecule object to calculate the energy of.
        fragment_indicies   - List of indicies of fragments to include in the calculation.
        model               - The model to use for this calculation, should be specified as method/basis.
        cp                  - Whether to use counterpoise correction, True means that fragments not included in
                    fragment_indicies will be added as ghost atoms to the calculation.
        settings            - SettingsReader object with all relevent settings.

    Returns:
        The calculated energy of the subset of fragments in this molecule.
    """

    method, basis = model.split("/")

    # Check if the user has qchem isntalled
    try:
        system.call("which", "qchem")
    except CommandExecutionError:
        raise LibraryNotAvailableError("qchem")

    # initialize qchem input string
    qchem_input = "";
    
    # molecule format
    qchem_input += "$molecule\n"

    # charge and spin multiplicity
    qchem_input += "{} {}\n".format(molecule.get_charge(), molecule.get_spin_multiplicity())

    # atoms in the molecule
    # might need to add whitespace before each line?
    qchem_input += molecule.to_xyz(fragment_indicies, cp) + "\n"

    qchem_input += "$end\n"

    # Q-chem settings
    qchem_input += "$rem\n"

    qchem_input += "jobtype " + "sp" + "\n"
    qchem_input += "method " + model.split('/')[0] + "\n"
    qchem_input += "basis " + model.split('/')[1] + "\n"

    try:
        qchem_input += "ecp " + settings.get("qchem", "ecp") + "\n"
    except ConfigMissingSectionError, ConfigMissingPropertyError:
        pass

    qchem_input += "$end"

    # file to write qchem input in
    qchem_in_path = files.get_energy_log_path(settings.get("files", "log_path"), molecule, method, basis, cp, "in")
    
    with open(qchem_in_path, "w") as qchem_in_file:
        f.write(qchem_input)
    
    # file to write qchem output in
    qchem_out_path = files.get_energy_log_path(settings.get("files", "log_path"), molecule, method, basis, cp, "out")

    # get number of threads
    num_threads = settings.getint("qchem", "num_threads")

    # perform system call to run qchem
    try:
        system.call("qchem", "-nt", num_threads, log_file_in, log_file_out)
    except CommandExecutionError as e:
        raise LibraryCallError("qchem", "energy calculation", "process returned non-zero exit code") from e

    # find the energy inside the qchem file output
    with open(qchem_out_path) as qchem_out_file:
        for line in qchem_out_file:
            if line.find("Total energy in the final basis set = ") != -1:
                return float(line[line.find("Total energy in the final basis set = ") + 39:])

    # if no line with "Total energy in the final basis set = " is found, raise an exception
    raise LibraryCallError("qchem", "energy calculation", "process returned file of incorrect format")

'''
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
'''
