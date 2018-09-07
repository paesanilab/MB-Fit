# qcalc.py
#
# Calculator that uses accepts calls from nmcgen and calls upon the requested quantum chemistry code (e.g. psi4, qchem, etc) to carry out the specified calculation.

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../qm_mb_energy_calculator/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import io, subprocess
import molecule_parser
from molecule import Molecule
from exceptions import ConfigMissingSectionError, ConfigMissingPropertyError
try:
    import psi4
except ImportError:
    has_psi4 = False
else:
    has_psi4 = True

def init(config, log_name):
    verify_program(config)
    
    if config['config_generator']['code'] == 'psi4':
        psi4_helper.init(config, log_name)

def optimize(settings, molecule, method, basis):
    
    if settings.get("config_generator", "code") == "psi4":
        return optimize_psi4(settings, molecule, method, basis)
    
    elif settings.get("config_generator", "code") == 'qchem':
        return optimize_qchem(settings, molecule, method, basis)

def optimize_psi4(settings, molecule, method, basis):

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
        raise

    energy = psi4.optimize("{}/{}".format(method, basis), molecule=psi4_mol)

    return Molecule().read_psi4_string(psi4_mol.save_string_xyz(), "molecule"), energy

def optimize_qchem(settings, molecule, method, basis):

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

    subprocess.run("qchem -nt {} {} {}".format(settings.getint("qchem", "num_threads", 1), qchem_input_path, qchem_output_path),
                   stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                   shell=True, check=True)

    found = False
    # found is set to true when the keyword 'Final energy is' is found. If found is true, it then looks for the keyword 'ATOM' to look for the optimized geometry
    qchem_output_string = ""
    with open(qchem_output_path, "r") as qchem_output_file:
        for line in qchem_output_file:
            if found:
                if "ATOM" in line:
                    for atom_index in range(molecule.get_num_atoms()):
                        qchem_output_string += " ".join(qchem_output_file.readline().split()[1:]) + "\n"
                    return Molecule().read_psi4_string("{} {}\n{}".format(molecule.get_charge(), molecule.get_spin_multiplicity(), qchem_output_string), "molecule"), energy
            elif "Final energy is" in line:
                energy = float(line.split()[3])
                found = True


    

def frequencies(settings, molecule, method, basis):

    if settings.get("config_generator", "code") == "psi4":
        return frequencies_psi4(settings, molecule, method, basis)
    
    elif settings.get("config_generator", "code") == 'qchem':
        return frequencies_qchem(settings, molecule, method, basis)

def frequencies_psi4(settings, molecule, method, basis):
    print("Using psi4 to determine normal modes.")

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
        raise

    total_energy, wavefunction = psi4.frequency("{}/{}".format(method, basis), molecule=psi4_mol, return_wfn=True)

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
    print("Using qchem to determine normal modes.")

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

    subprocess.run("qchem -nt {} {} {}".format(settings.getint("qchem", "num_threads", 1), qchem_input_path, qchem_output_path), 
                   stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                   shell=True, check=True)

    num_atoms = molecule.get_num_atoms()

    with open(qchem_output_path, "r") as qchem_output_file:
        frequencies = []
        red_masses = []
        normal_modes = []
        while(True):
            line = qchem_output_file.readline()

            if line == "":
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
