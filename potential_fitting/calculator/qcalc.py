# external package imports
import subprocess, os

# absolute module imports
from potential_fitting.utils import files, system
from potential_fitting.molecule import Molecule
from potential_fitting.exceptions import (LibraryNotAvailableError, NoSuchLibraryError, ConfigMissingSectionError,
        ConfigMissingPropertyError, CommandNotFoundError, CommandExecutionError, LibraryCallError)

try:
    import psi4
    from psi4.driver.qcdb.exceptions import QcdbException
except ImportError:
    has_psi4 = False
else:
    has_psi4 = True

def optimize(settings, molecule, method, basis):
    """
    Optimizes the geometry of the given molecule with the given method and basis.

    The input molecule will be unchanged.

    Args:
        settings            - SettingsReader object with relevant settings.
        molecule            - The Molecule object to optimize.
        method              - The method to optimize with.
        basis               - The basis to optimize with.

    Returns:
        A 2-tuple, containing a new Molecule object of the optimized geometry and the energy of the optimized geometry.
    """
    
    library = settings.get("config_generator", "code")

    # check if the user requested a psi4 calculation
    if library == "psi4":
        return optimize_psi4(settings, molecule, method, basis)
    
    # check if the user requested a qchem calculation
    elif library == "qchem":
        return optimize_qchem(settings, molecule, method, basis)

    # raise an exception if the user requested an unsupported library
    raise NoSuchLibraryError(library)

def optimize_psi4(settings, molecule, method, basis):
    """
    Optimizes the geometry of the given molecule using psi4 with the given method and basis.

    The input molecule will be unchanged.

    Args:
        settings            - SettingsReader object with relevant settings.
        molecule            - The Molecule object to optimize.
        method              - The method to optimize with.
        basis               - The basis to optimize with.

    Returns:
        A 2-tuple, containing a new Molecule object of the optimized geometry and the energy of the optimized geometry.
    """

    # raise an exception if psi4 is not installed
    if not has_psi4:
        raise LibraryNotAvailableError("psi4")

    print("Beginning geometry optimization using psi4 of {} with {}/{}.".format(molecule.get_name(), method, basis))

    log_path = files.get_optimization_log_path(settings.get("files", "log_path"), molecule, method, basis, "log")

    # set psi4's log path, memory, and number of threads to use.
    psi4.core.set_output_file(log_path, False)
    psi4.set_memory(settings.get("psi4", "memory"))
    psi4.set_num_threads(settings.getint("psi4", "num_threads"))

    # construct an input string for psi4.geometry()
    psi4_string = "{}\n{} {}".format(molecule.to_xyz(), molecule.get_charge(), molecule.get_spin_multiplicity())

    # construct a psi4 molecule from the input string
    try:
        psi4_mol = psi4.geometry(psi4_string)
    except RuntimeError as e:
        raise LibraryCallError("psi4", "geometry", str(e))

    # use psi4 to optimize the geometry
    try:
        energy = psi4.optimize("{}/{}".format(method, basis), molecule=psi4_mol)
    except QcdbException as e:
        raise LibraryCallError("psi4", "optimize", str(e))

    print("Completed geometry optimization.")

    return Molecule().read_psi4_string(psi4_mol.save_string_xyz()), energy

def optimize_qchem(settings, molecule, method, basis):
    """
    Optimizes the geometry of the given molecule using qchem with the given method and basis.

    The input molecule will be unchanged.

    Args:
        settings            - SettingsReader object with relevant settings.
        molecule            - The Molecule object to optimize.
        method              - The method to optimize with.
        basis               - The basis to optimize with.

    Returns:
        A 2-tuple, containing a new Molecule object of the optimized geometry and the energy of the optimized geometry.
    """

    # Check if the user has qchem isntalled
    try:
        system.call("which", "qchem")
    except CommandExecutionError:
        raise LibraryNotAvailableError("qchem")

    print("Beginning geometry optimization using qchem of {} with {}/{}.".format(molecule.get_name(), method, basis))

    qchem_in_path = files.get_optimization_log_path(settings.get("files", "log_path"), molecule, method, basis, "in")

    # construct the qchem input file
    with open(qchem_in_path, "w") as qchem_in_file:

        # tells qchem that the molecule starts here
        qchem_in_file.write("$molecule\n")
        # this line contains charge and spin multiplicity
        qchem_in_file.write("{} {}\n".format(molecule.get_charge(), molecule.get_spin_multiplicity()))
        # write the molecule's xyz format
        qchem_in_file.write(molecule.to_xyz())
        # tells qchem that the molecule ends here
        qchem_in_file.write("\n$end\n")

        qchem_in_file.write("\n")

        # tells qchem that configuration starts here
        qchem_in_file.write("$rem\n")
        # tells qchem this is an optimization job
        qchem_in_file.write("jobtype opt\n")

        # tells qchem what method and basis to use
        qchem_in_file.write("method " + method + "\n")
        qchem_in_file.write("basis " + basis + "\n")
        qchem_in_file.write("GEOM_OPT_MAX_CYCLES 200\n")

        try:
            qchem_in_file.write("ecp {}\n".format(settings.get("config_generator", "ecp")))
        except (ConfigMissingSectionError, ConfigMissingPropertyError):
            pass

        # tells qchem that configuration ends here
        qchem_in_file.write("$end\n")

    qchem_out_path = files.get_optimization_log_path(settings.get("files", "log_path"), molecule, method, basis, "out")

    # make the qchem system call
    try:
        system.call("qchem", "-nt", str(settings.getint("qchem", "num_threads", 1)), qchem_in_path, qchem_out_path)
    except CommandExecutionError as e:
        raise LibraryCallError("qchem", "optimize", "process returned non-zero exit code") from e
    
    # parse the molecule's geometry out of the output file
    found = False
    # found is set to true when the keyword 'Final energy is' is found. If found is true, it then looks for the keyword
    # 'ATOM' to look for the optimized geometry
    qchem_out_string = ""
    with open(qchem_out_path, "r") as qchem_out_file:
        for line in qchem_out_file:
            if found:
                if "ATOM" in line:
                    for atom_index in range(molecule.get_num_atoms()):
                        qchem_out_string += " ".join(qchem_out_file.readline().split()[1:]) + "\n"
                    print("Completed geometry optimization.")
                    return Molecule().read_psi4_string("{} {}\n{}".format(molecule.get_charge(),
                            molecule.get_spin_multiplicity(), qchem_out_string)), energy
            elif "Final energy is" in line:
                energy = float(line.split()[3])
                found = True

    # if we didn't find a geometry to parse, raise an exception.
    raise LibraryCallError("qchem", "optimze", "process returned file of incorrect format")

def frequencies(settings, molecule, method, basis):
    """
    Performs a frequency calculation on the molecule with the given method and basis.

    Args:
        settings            - SettingsReader object with relevant settings.
        molecule            - The Molecule object to perform a frequency calculation on.
        method              - The method to use for the frequency calculation.
        basis               - The basis to use for the frequency calculation.

    Returns:
        A 3-tuple containing the normal modes, frequencies, and reduced masses of the molecule.
    """

    if settings.get("config_generator", "code") == "psi4":
        return frequencies_psi4(settings, molecule, method, basis)
    
    elif settings.get("config_generator", "code") == 'qchem':
        return frequencies_qchem(settings, molecule, method, basis)

    # raise an exception if the user requested an unsupported library
    raise NoSuchLibraryError(library)

def frequencies_psi4(settings, molecule, method, basis):
    """
    Performs a frequency calculation on the molecule using psi4 with the given method and basis.

    Args:
        settings            - SettingsReader object with relevant settings.
        molecule            - The Molecule object to perform a frequency calculation on.
        method              - The method to use for the frequency calculation.
        basis               - The basis to use for the frequency calculation.

    Returns:
        A 3-tuple containing the normal modes, frequencies, and reduced masses of the molecule.
    """

    print("Beginning normal modes calculation using psi4 of {} with {}/{}.".format(molecule.get_name(), method, basis))

    log_path = files.get_frequencies_log_path(settings.get("files", "log_path"), molecule, method, basis, "log")

    # set psi4's log path, memory, and number of threads to use.
    psi4.core.set_output_file(log_path, False)
    psi4.set_memory(settings.get("psi4", "memory"))
    psi4.set_num_threads(settings.getint("psi4", "num_threads"))

    # construct an input string for psi4.geometry()
    psi4_string = "{}\n{} {}".format(molecule.to_xyz(), molecule.get_charge(), molecule.get_spin_multiplicity())

    # construct a psi4 molecule from the input string
    try:
        psi4_mol = psi4.geometry(psi4_string)
    except RuntimeError as e:
        raise LibraryCallError("psi4", "geometry", str(e))

    # use psi4 to perform a frequency calculation
    try:
        total_energy, wavefunction = psi4.frequency("{}/{}".format(method, basis), molecule=psi4_mol, return_wfn=True)
    except QcdbException as e:
        raise LibraryCallError("psi4", "frequency", str(e))

    # retrieve the normal modes, frequencies, and reduced masses from the psi4 output object
    vib_info = psi4.qcdb.vib.filter_nonvib(wavefunction.frequency_analysis)

    frequencies = wavefunction.frequencies().to_array()
    red_masses = vib_info["mu"].data
    
    normal_modes_raw = vib_info['x'].data
    normal_modes = []
    
    num_modes = len(frequencies)
    num_atoms = molecule.get_num_atoms()

    for normal_mode_index in range(num_modes):
        normal_mode = [[0, 0, 0] for atom_index in range(num_atoms)]
        
        for atom_index in range(num_atoms):
            normal_mode[atom_index][0] = float(normal_modes_raw[3 * atom_index + 0][normal_mode_index])
            normal_mode[atom_index][1] = float(normal_modes_raw[3 * atom_index + 1][normal_mode_index])
            normal_mode[atom_index][2] = float(normal_modes_raw[3 * atom_index + 2][normal_mode_index])
        
        normal_modes.append(normal_mode)

    print("Normal mode/frequency analysis complete. {} normal modes found.".format(num_modes))

    num_neg_freqs = 0
    for frequency in frequencies:
        if frequency < 0:
            num_neg_freqs += 1

    if num_neg_freqs == 1:
        print("Single negative frequency detected. This most likely means the given geometry is a transition state.")
        
    elif num_neg_freqs > 1:
        print("Multiple ({}) negative frequencies detected. Proceed with caution.".format(num_neg_freqs))
    
    return normal_modes, frequencies, red_masses


def frequencies_qchem(settings, molecule, method, basis):
    """
    Performs a frequency calculation on the molecule using qchem with the given method and basis.

    Args:
        settings            - SettingsReader object with relevant settings.
        molecule            - The Molecule object to perform a frequency calculation on.
        method              - The method to use for the frequency calculation.
        basis               - The basis to use for the frequency calculation.

    Returns:
        A 3-tuple containing the normal modes, frequencies, and reduced masses of the molecule.
    """

    # Check if the user has qchem isntalled
    try:
        system.call("which", "qchem")
    except CommandExecutionError:
        raise LibraryNotAvailableError("qchem")

    print("Beginning normal modes calculation using qchem of {} with {}/{}.".format(molecule.get_name(), method,
            basis))

    qchem_in_path = files.get_frequencies_log_path(settings.get("files", "log_path"), molecule, method, basis, "in")

    with open(qchem_in_path, "w") as qchem_in_file:

        # tells qchem that the molecule starts here
        qchem_in_file.write("$molecule\n")
        # this line contains charge and spin multiplicity
        qchem_in_file.write("{} {}\n".format(molecule.get_charge(), molecule.get_spin_multiplicity()))
        # write the molecule's xyz format
        qchem_in_file.write(molecule.to_xyz())
        # tells qchem that the molecule ends here
        qchem_in_file.write("$end\n")

        qchem_in_file.write("\n")

        # tells qchem that configuration starts here
        qchem_in_file.write("$rem\n")
        # tells qchem this is an optimization job
        qchem_in_file.write("jobtype freq\n")

        # tells qchem what method and basis to use
        qchem_in_file.write("method " + method + "\n")
        qchem_in_file.write("basis " + basis + "\n")

        try:
            qchem_in_file.write("ecp " + settings.get("config_generator", "ecp") + "\n")
        except (ConfigMissingSectionError, ConfigMissingPropertyError):
            pass

        # tells qchem that configuration ends here
        qchem_in_file.write("$end\n")

    qchem_out_path = files.get_frequencies_log_path(settings.get("files", "log_path"), molecule, method, basis, "out")

    try:
        system.call("qchem", "-nt", str(settings.getint("qchem", "num_threads", 1)), qchem_in_path, qchem_out_path)
    except CommandExecutionError:
        raise LibraryCallError("qchem", "frequency", "process returned non-zero exit code")

    num_atoms = molecule.get_num_atoms()

    # now parse the qchem output file to read the normal modes, frequencies, and reduced masses
    with open(qchem_out_path, "r") as qchem_out_file:
        frequencies = []
        red_masses = []
        normal_modes = []
        while(True):
            line = qchem_out_file.readline()

            if line == "":
                if frequencies == []:
                    raise LibraryCallError("qchem", "frequency", "process returned file of incorrect format")
                break

            if "Mode:" in line:
                num_modes = len(line.split()) - 1

                # read frequencies from file and cast each to a float
                frequencies += [float(freq) for freq in qchem_out_file.readline().split()[1:]]

                # skip the Force Constant line
                qchem_out_file.readline()

                # read red masses from file and cast each to float
                red_masses += [float(freq) for freq in qchem_out_file.readline().split()[2:]]

                # skip IR Active, IR Intens, Raman Active and labels line
                for i in range(4):
                    qchem_out_file.readline()

                # initialize 3d list for normal modes
                modes = [[[0, 0, 0] for i in range(num_atoms)] for k in range(num_modes)]

                for atom in range(num_atoms):

                    # read the coordinates for the next atom for the next num_modes modes
                    modes_of_atom = [float(i) for i in qchem_out_file.readline().split()[1:]]

                    # loop over each mode
                    for index in range(num_modes):

                        # store coordinates in corresponding entries of modes list
                        modes[index][atom][0] = modes_of_atom[index * 3]
                        modes[index][atom][1] = modes_of_atom[index * 3 + 1]
                        modes[index][atom][2] = modes_of_atom[index * 3 + 2]

                normal_modes += modes

    print("Normal mode/frequency analysis complete. {} normal modes found".format(len(normal_modes)))

    num_neg_freqs = 0
    for frequency in frequencies:
        if frequency < 0:
            num_neg_freqs += 1

    if num_neg_freqs == 1:
        print("Single negative frequency detected. This most likely means the given geometry is a transition state.")
    
    elif num_neg_freqs > 1:
        print("Multiple ({}) negative frequencies detected. Proceed with caution.".format(num_neg_freqs))
    
    return normal_modes, frequencies, red_masses
