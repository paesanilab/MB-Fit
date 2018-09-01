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

    subprocess.run("qchem -nt {} {} {}".format(settings.getint("qchem", "num_threads"), qchem_input_path, qchem_output_path),
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


    

def frequencies(molecule, filenames, config):

    verify_program(config)

    if config['config_generator']['code'] == 'psi4':   
        psi4_molecule = psi4_helper.psi4_mol(molecule, config['molecule']['charges'], config['molecule']['spins'])
        return psi4_helper.frequencies(psi4_molecule, config)

    elif config['config_generator']['code'] == 'qchem':  
        return qchem_helper.frequencies(molecule, filenames, config)

