# qcalc.py
#
# Calculator that uses accepts calls from nmcgen and calls upon the requested quantum chemistry code (e.g. psi4, qchem, etc) to carry out the specified calculation.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../qm_mb_energy_calculator/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import io, subprocess
import molecule_parser
from molecule import Molecule
from exceptions import LibraryNotAvailableError, NoSuchLibraryError, ConfigMissingSectionError, ConfigMissingPropertyError
from psi4.driver.qcdb.exceptions import QcdbException
try:
    import psi4
except ImportError:
    has_psi4 = False
else:
    has_psi4 = True

def optimize(settings, molecule, method, basis):
    
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

    # raise an exception if psi4 is not installed
    if not has_psi4:
        raise LibraryNotAvailableError("psi4")

    print("Beginning geometry optimization using psi4 of {} with {}/{}.".format(molecule.get_name(), method, basis))

    log_path = settings.get("files", "log_path") + "/optimizations/{}/{}/{}.out".format(method, basis, molecule.get_SHA1()[-8:])

    if not os.path.exists(os.path.dirname(log_path)):
        os.makedirs(os.path.dirname(log_path))

    psi4.core.set_output_file(log_path, False)
    psi4.set_memory(settings.get("psi4", "memory"))
    psi4.set_num_threads(settings.getint("psi4", "num_threads"))

    psi4_string = "{}\n{} {}".format(molecule.to_xyz(), molecule.get_charge(), molecule.get_spin_multiplicity())

    try:
        psi4_mol = psi4.geometry(psi4_string)
    except RuntimeError as e:
        raise LibraryCallError("psi4", "geometry", str(e))

    try:
        energy = psi4.optimize("{}/{}".format(method, basis), molecule=psi4_mol)
    except QcdbException as e:
        raise LibraryCallerror("psi4", "optimize", str(e))

    print("Completed geometry optimization.")

    return Molecule().read_psi4_string(psi4_mol.save_string_xyz()), energy

def optimize_qchem(settings, molecule, method, basis):

    print("Beginning geometry optimization using qchem of {} with {}/{}/".format(molecule.get_name(), method, basis))

    qchem_input_path = settings.get("files", "log_path") + "/optimizations/{}/{}/{}.in".format(method, basis, molecule.get_SHA1()[-8:])

    if not os.path.exists(os.path.dirname(qchem_input_path)):
        os.makedirs(os.path.dirname(qchem_input_path))

    with open(qchem_input_path, "w") as qchem_input_file:

        # tells qchem that the molecule starts here
        qchem_input_file.write("$molecule\n")
        # this line contains charge and spin multiplicity
        qchem_input_file.write("{} {}\n".format(molecule.get_charge(), molecule.get_spin_multiplicity()))
        # write the molecule's xyz format
        qchem_input_file.write(molecule.to_xyz())
        # tells qchem that the molecule ends here
        qchem_input_file.write("$end\n")

        qchem_input_file.write("\n")

        # tells qchem that configuration starts here
        qchem_input_file.write("$rem\n")
        # tells qchem this is an optimization job
        qchem_input_file.write("jobtype opt\n")

        # tells qchem what method and basis to use
        qchem_input_file.write("method " + method + "\n")
        qchem_input_file.write("basis " + basis + "\n")

        try:
            qchem_input_file.write("ecp " + settings.get("config_generator", "ecp") + "\n")
        except (ConfigMissingSectionError, ConfigMissingPropertyError):
            pass

        # tells qchem that configuration ends here
        qchem_input_file.write("$end\n")
    qchem_output_path = settings.get("files", "log_path") + "/optimizations/{}/{}/{}.out".format(method, basis, molecule.get_SHA1()[-8:])

    if not os.path.exists(os.path.dirname(qchem_output_path)):
        os.makedirs(os.path.dirname(qchem_output_path))

    qchem_call_string = "qchem -nt {} {} {}".format(settings.getint("qchem", "num_threads", 1), qchem_input_path, qchem_output_path)
    if subprocess.run(qchem_call_string,
                   stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                   shell=True, check=True).returncode != 0:
        raise LibraryCallError("qchem", "optimize", "process returned non-zero exit code")
        

    found = False
    # found is set to true when the keyword 'Final energy is' is found. If found is true, it then looks for the keyword 'ATOM' to look for the optimized geometry
    qchem_output_string = ""
    with open(qchem_output_path, "r") as qchem_output_file:
        for line in qchem_output_file:
            if found:
                if "ATOM" in line:
                    for atom_index in range(molecule.get_num_atoms()):
                        qchem_output_string += " ".join(qchem_output_file.readline().split()[1:]) + "\n"
                    print("Completed geometry optimization.")
                    return Molecule().read_psi4_string("{} {}\n{}".format(molecule.get_charge(), molecule.get_spin_multiplicity(), qchem_output_string)), energy
            elif "Final energy is" in line:
                energy = float(line.split()[3])
                found = True

    raise LibraryCallError("qchem", "optimze", "process returned file of incorrect format")

def frequencies(settings, molecule, method, basis):

    if settings.get("config_generator", "code") == "psi4":
        return frequencies_psi4(settings, molecule, method, basis)
    
    elif settings.get("config_generator", "code") == 'qchem':
        return frequencies_qchem(settings, molecule, method, basis)

    # raise an exception if the user requested an unsupported library
    raise NoSuchLibraryError(library)

def frequencies_psi4(settings, molecule, method, basis):
    print("Beginning normal modes calculation using psi4 of {} with {}/{}/".format(molecule.get_name(), method, basis))

    log_path = settings.get("files", "log_path") + "/normal_modes/{}/{}/{}.out".format(method, basis, molecule.get_SHA1()[-8:])

    if not os.path.exists(os.path.dirname(log_path)):
        os.makedirs(os.path.dirname(log_path))

    psi4.core.set_output_file(log_path, False)

    psi4.set_memory(settings.get("psi4", "memory"))
    psi4.set_num_threads(settings.getint("psi4", "num_threads"))

    psi4_string = "{}\n{} {}".format(molecule.to_xyz(), molecule.get_charge(), molecule.get_spin_multiplicity())

    try:
        psi4_mol = psi4.geometry(psi4_string)
    except RuntimeError as e:
        raise LibraryCallError("psi4", "geometry", str(e))

    try:
        total_energy, wavefunction = psi4.frequency("{}/{}".format(method, basis), molecule=psi4_mol, return_wfn=True)
    except QcdbException as e:
        raise LibraryCallError("psi4", "frequency", str(e))

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

    print("Normal mode/frequency analysis complete. {} normal modes found".format(num_modes))
    
    return normal_modes, frequencies, red_masses


def frequencies_qchem(settings, molecule, method, basis):
    print("Beginning normal modes calculation using qchem of {} with {}/{}/".format(molecule.get_name(), method, basis))

    qchem_input_path = settings.get("files", "log_path") + "/normal_modes/{}/{}/{}.in".format(method, basis, molecule.get_SHA1()[-8:])

    if not os.path.exists(os.path.dirname(qchem_input_path)):
        os.makedirs(os.path.dirname(qchem_input_path))

    with open(qchem_input_path, "w") as qchem_input_file:

        # tells qchem that the molecule starts here
        qchem_input_file.write("$molecule\n")
        # this line contains charge and spin multiplicity
        qchem_input_file.write("{} {}\n".format(molecule.get_charge(), molecule.get_spin_multiplicity()))
        # write the molecule's xyz format
        qchem_input_file.write(molecule.to_xyz())
        # tells qchem that the molecule ends here
        qchem_input_file.write("$end\n")

        qchem_input_file.write("\n")

        # tells qchem that configuration starts here
        qchem_input_file.write("$rem\n")
        # tells qchem this is an optimization job
        qchem_input_file.write("jobtype freq\n")

        # tells qchem what method and basis to use
        qchem_input_file.write("method " + method + "\n")
        qchem_input_file.write("basis " + basis + "\n")

        try:
            qchem_input_file.write("ecp " + settings.get("config_generator", "ecp") + "\n")
        except (ConfigMissingSectionError, ConfigMissingPropertyError):
            pass

        # tells qchem that configuration ends here
        qchem_input_file.write("$end\n")
    qchem_output_path = settings.get("files", "log_path") + "/normal_modes/{}/{}/{}.out".format(method, basis, molecule.get_SHA1()[-8:])

    if not os.path.exists(os.path.dirname(qchem_output_path)):
        os.makedirs(os.path.dirname(qchem_output_path))

    qchem_call_string = "qchem -nt {} {} {}".format(settings.getint("qchem", "num_threads", 1), qchem_input_path, qchem_output_path)
    if subprocess.run(qchem_call_string, 
                   stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                   shell=True, check=True).returncode != 0:
        raise LibraryCallError("qchem", "frequency", "process returned non-zero exit code")

    num_atoms = molecule.get_num_atoms()

    with open(qchem_output_path, "r") as qchem_output_file:
        frequencies = []
        red_masses = []
        normal_modes = []
        while(True):
            line = qchem_output_file.readline()

            if line == "":
                if frequencies == []:
                    raise LibraryCallError("qchem", "frequency", "process returned file of incorrect format")
                break

            if "Mode:" in line:
                num_modes = len(line.split()) - 1

                # read frequencies from file and cast each to a float
                frequencies += [float(freq) for freq in qchem_output_file.readline().split()[1:]]

                # skip the Force Constant line
                qchem_output_file.readline()

                # read red masses from file and cast each to float
                red_masses += [float(freq) for freq in qchem_output_file.readline().split()[2:]]

                # skip IR Active, IR Intens, Raman Active and labels line
                for i in range(4):
                    qchem_output_file.readline()

                # initialize 3d list for normal modes
                modes = [[[0, 0, 0] for i in range(num_atoms)] for k in range(num_modes)]

                for atom in range(num_atoms):

                    # read the coordinates for the next atom for the next num_modes modes
                    modes_of_atom = [float(i) for i in qchem_output_file.readline().split()[1:]]

                    # loop over each mode
                    for index in range(num_modes):

                        # store coordinates in corresponding entries of modes list
                        modes[index][atom][0] = modes_of_atom[index * 3]
                        modes[index][atom][1] = modes_of_atom[index * 3 + 1]
                        modes[index][atom][2] = modes_of_atom[index * 3 + 2]

                normal_modes += modes

    print("Normal mode/frequency analysis complete. {} normal modes found".format(len(normal_modes)))

    return normal_modes, frequencies, red_masses
