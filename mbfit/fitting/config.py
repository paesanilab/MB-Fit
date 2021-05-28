import os,sys
import math, configparser, itertools
from collections import OrderedDict

from mbfit.utils import constants, SettingsReader, files, system
from mbfit.exceptions import InvalidValueError, InconsistentValueError, LibraryCallError
from mbfit.molecule import Molecule, Fragment, xyz_to_molecules
from mbfit.polynomials import MoleculeInParser

qchem_template = "qchem_template"

# TODO This function should be removed and use the fancier fucntion that Ethan has

# to get the list of atoms
def get_atom_types(fragment):
    atom_list = []
    was_digit = False
    current_text = ""
    for character in fragment:
        if not character.isdigit():
            if was_digit:
                atom_list.append(int(current_text))
                current_text = character
                was_digit = False
            else:
                current_text += character
        else:
            if was_digit:
                current_text += character
            else:
                atom_list.append(current_text)
                was_digit = True
                current_text = character

    # At this point, only the last number is missing
    atom_list.append(int(current_text))
    return atom_list

def calculate_c6_for_config(dimer_settings, settings_list, geo_list, fragment_list, distance_between, use_published_polarizabilities,
                            method,
                            basis, use_cm5):
    # Define the names
    name1 = settings_list[0].get("molecule","names")
    name2 = settings_list[1].get("molecule","names")

    print("Running C6 calculation for the dimer {}_{} with {}/{}.".format(name1, name2, method, basis))


    qchem_id_name = "get_config_qchem_{}_{}_{}_{}".format(name1, name2, method, basis)

    qchem_in_path = files.init_file(os.path.join(dimer_settings.get("files", "log_path"), qchem_id_name + ".in"), files.OverwriteMethod.OVERWRITE)
    qchem_out_path = files.init_file(os.path.join(dimer_settings.get("files", "log_path"), qchem_id_name + ".out"), files.OverwriteMethod.OVERWRITE)
    qchem_log_path = files.init_file(os.path.join(dimer_settings.get("files", "log_path"), qchem_id_name + ".log"), files.OverwriteMethod.OVERWRITE)

    # read geometry into a molecule object
    molecule1 = xyz_to_molecules(geo_list[0], settings=settings_list[0])[0]
    molecule2 = xyz_to_molecules(geo_list[1], settings=settings_list[1])[0]

    # move molecule1 to its standard orientation
    molecule1.move_to_center_of_mass()
    molecule1.rotate_on_principal_axes()

    # Same with mol 2, but also move it away
    molecule2.move_to_center_of_mass()
    molecule2.rotate_on_principal_axes()
    molecule2.translate(distance_between, 0, 0)

    # Get the atomic symbols list
    atomic_symbols = []
    for atom in molecule1.get_atoms():
        atomic_symbols.append(atom.get_name())
    for atom in molecule2.get_atoms():
        atomic_symbols.append(atom.get_name())

    if os.path.exists(qchem_out_path):
        print("Skipping dimer {}_{}. This calculation has already been performed.".format(name1,name2))
    else:

        with open(qchem_in_path, "w") as qchem_in:
            # tells qchem that the molecule has started
            qchem_in.write("$molecule\n")
            # charge and spin line
            qchem_in.write("{} {}\n".format(sum([int(charge) for charge in dimer_settings.get("molecule", "charges").split(",")]),
                                            1 + sum(
                                                [int(spin) - 1 for spin in dimer_settings.get("molecule", "spins").split(",")])))
            # copy this molecule to the qchem input
            qchem_in.write(molecule1.to_xyz())
            qchem_in.write("\n")
            qchem_in.write(molecule2.to_xyz())
            qchem_in.write("\n")

            # tells qchem that the molecule has ended
            qchem_in.write("$end\n")

            qchem_in.write("$rem\n")
            qchem_in.write("method " + method + "\n")
            qchem_in.write("basis " + basis + "\n")

            # read the qchem template and append it to the qchem in
            qchem_template2 = qchem_template 
            if use_cm5 == False:
                qchem_template2 += "_chelpg"
            qchem_template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), qchem_template2)
            with open(qchem_template_path, "r") as template:
                for line in template:
                    qchem_in.write(line)

            num_threads = dimer_settings.getint("qchem", "num_threads", 1)

        with open(qchem_log_path, "w") as qchem_log:
            system.call("qchem", "-nt", str(num_threads), qchem_in_path, qchem_out_path, out_file=qchem_log)

    print("Parsing qchem output {}.".format(qchem_out_path))

    c6 = get_c6_from_qchem_output(qchem_out_path, fragment_list, atomic_symbols, use_published_polarizabilities)

    return c6

def get_c6_from_qchem_output(qchem_out_path, fragments, atomic_symbols, use_published_polarizabilities):
    # parse the output file
    with open(qchem_out_path, "r") as qchem_out:
        # read lines until we read the line before where c6 constants are specified
        while True:
            line = qchem_out.readline()
            if "SCF failed to converge" in line:
                raise LibraryCallError("Qchem", "c6 calculation", "SCF failed to converge".format(qchem_out_path), log_path = qchem_out_path)
            try:
                line.index("D2X   (i,j)")
                break
            except ValueError:
                pass

        # each element of this dictionary stores all the c6 constants for a single fragment, with the last dictionary storing inter-molecular c6 constants
        # each dictionary is a dictionary of arrays, with the mapping being for example "AB" -> [array of all AB c6] this allows us to average all the AB c6
        # constants later to get a more accurate c6. Right now, inter-molecular c6 are treated differently, and are averaged seperately from intra-molecular c6
        c6_constant_lists = [OrderedDict() for x in range(len(fragments) + 1)]

        # keeps track of which fragment atom a is in (ie 0->A3B2C1, 1->D2 for A3B2C1_D2)
        fragment_index_a = 0

        # keeps track of which letter corresponds to atom a (ie 0->A, 1->B, 2->C for A3B2C1)
        atom_index_a = 0

        # keeps track of which atom within a letter corresponds to atom a (ie 0, then 1, then 2 for A3)
        atom_a = 0

        while (True):
            # check if atom_a is out of bounds of the number of atoms for the current a atom type
            if atom_a >= int(fragments[fragment_index_a][atom_index_a * 2 + 1]):

                # move to the next atom a type
                atom_index_a += 1
                atom_a = 0

                # check if the atom type a is out of bounds for the current fragment a
                if atom_index_a >= len(fragments[fragment_index_a][::2]):

                    # move to the next fragment a
                    fragment_index_a += 1
                    atom_index_a = 0

                    # check if the current fragment a is out of bounds, in which case all c6 constants have been parsed
                    if fragment_index_a >= len(fragments):
                        break

            # keeps track of which fragment b is in
            fragment_index_b = fragment_index_a

            # keeps track of which letter corresponds to atom b
            atom_index_b = atom_index_a

            # keeps track of which atom within a letter corresponds to atom b, starts 1 atom after atom a so each combination of atom is iterated over once.
            atom_b = atom_a + 1

            while (True):
                # check if atom_b is out of bounds of the number of atoms for the current b atom type
                if atom_b >= int(fragments[fragment_index_b][atom_index_b * 2 + 1]):

                    # move to the next atom b type
                    atom_index_b += 1
                    atom_b = 0

                    # check if the atom type b is out of bounds for the current fragment b
                    if atom_index_b >= len(fragments[fragment_index_b][::2]):

                        # move to the next fragment b
                        fragment_index_b += 1
                        atom_index_b = 0

                        # check if the current fragment b is out of bounds, in which case all c6 constants corresponding to atom a have been parsed
                        if fragment_index_b >= len(fragments):
                            break

                # read a new line from the inpit
                line = qchem_out.readline()

                # parse the c6 value from the line
                c6 = float(line.split()[2])

                # convert c6 constant from Haretrees*(Bohr^6) to (kcal/mol)*A^6
                c6 *= constants.au_times_bohr6_to_kcal_times_ang6

                # get each atom type based on the current fragments and atom type indicies
                atom_type_a = fragments[fragment_index_a][::2][atom_index_a]
                atom_type_b = fragments[fragment_index_b][::2][atom_index_b]

                atoms_a_first = atom_type_a + atom_type_b
                atoms_b_first = atom_type_b + atom_type_a

                # check if this is an intrafragmental c6 constant
                if fragment_index_a == fragment_index_b:
                    # loop over all fragments and add this c6 constant to any dictionary of an equivelent fragment
                    for fragment_index in range(len(fragments)):
                        if fragments[fragment_index] == fragments[fragment_index_a]:
                            if atoms_a_first in c6_constant_lists[fragment_index]:
                                c6_constant_lists[fragment_index][atoms_a_first].append(c6)

                            elif atoms_b_first in c6_constant_lists[fragment_index]:
                                c6_constant_lists[fragment_index][atoms_b_first].append(c6)

                            else:
                                c6_constant_lists[fragment_index][atoms_a_first] = []
                                c6_constant_lists[fragment_index][atoms_a_first].append(c6)

                # otherwise this is an interfragmental c6 constant
                else:

                    # add this c6 constant to the last dictionary (the one for interfragmental c6 constants)
                    if atoms_a_first in c6_constant_lists[len(fragments)]:
                        c6_constant_lists[len(fragments)][atoms_a_first].append(c6)

                    elif atoms_b_first in c6_constant_lists[len(fragments)]:
                        c6_constant_lists[len(fragments)][atoms_b_first].append(c6)

                    else:
                        c6_constant_lists[len(fragments)][atoms_a_first] = []
                        c6_constant_lists[len(fragments)][atoms_a_first].append(c6)

                atom_b += 1

            atom_a += 1
        # now construct the c6 constants 2d array by averagining all equivelent c6 constants,

        c6_constants = []

        for c6_dictionary in c6_constant_lists:
            c6_consts = []

            for key in c6_dictionary:
                c6_consts.append(sum(c6_dictionary[key]) / len(c6_dictionary[key]))

            c6_constants.append(c6_consts)

    return c6_constants



def calculate_chg_pol_for_config(settings, geo, fragment, distance_between, use_published_polarizabilities,
                                 method,
                                 basis,use_cm5):
    name = settings.get("molecule","names")
    print("Running C6 calculation for the monomer {} with {}/{}.".format(name, method, basis))

    qchem_id_name = "get_config_qchem_{}_{}_{}".format(name, method, basis)

    qchem_out_path = files.init_file(os.path.join(settings.get("files", "log_path"), qchem_id_name + ".out"), files.OverwriteMethod.OVERWRITE)
    qchem_in_path = files.init_file(os.path.join(settings.get("files", "log_path"), qchem_id_name + ".in"), files.OverwriteMethod.OVERWRITE)
    qchem_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), qchem_id_name + ".log"), files.OverwriteMethod.OVERWRITE)

    atomic_symbols = []

    # read geometry into a molecule object
    molecule = xyz_to_molecules(geo, settings)[0]

    # move molecule1 to its standard orientation
    molecule.move_to_center_of_mass()
    molecule.rotate_on_principal_axes()

    # add molecule1's atoms to the list of atomic_symbols
    for atom in molecule.get_atoms():
        atomic_symbols.append(atom.get_name())

    if os.path.exists(qchem_out_path):
        print("Skipping monomer {}. This calculation has already been performed.".format(name))
    else:
        with open(qchem_in_path, "w") as qchem_in:
            # tells qchem that the molecule has started
            qchem_in.write("$molecule\n")
            # charge and spin line
            qchem_in.write("{} {}\n".format(sum([int(charge) for charge in settings.get("molecule", "charges").split(",")]),
                                            1 + sum(
                                                [int(spin) - 1 for spin in settings.get("molecule", "spins").split(",")])))
            # copy this molecule to the qchem input
            qchem_in.write(molecule.to_xyz())
            qchem_in.write("\n")

            # tells qchem that the molecule has ended
            qchem_in.write("$end\n")

            qchem_in.write("$rem\n")
            qchem_in.write("method " + method + "\n")
            qchem_in.write("basis " + basis + "\n")

            # read the qchem template and append it to the qchem in
            qchem_template2 = qchem_template
            if use_cm5 == False:
                qchem_template2 += "_chelpg"
            qchem_template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), qchem_template2)
            with open(qchem_template_path, "r") as template:
                for line in template:
                    qchem_in.write(line)

            num_threads = settings.getint("qchem", "num_threads", 1)

        print("Executing qchem calculation...")
        # perform qchem system call
        with open(qchem_log_path, "w") as qchem_log:
            system.call("qchem", "-nt", str(num_threads), qchem_in_path, qchem_out_path, out_file=qchem_log)

    print("Parsing qchem output {}.".format(qchem_out_path))

    charges, pols = get_chg_pol_from_qchem_output(qchem_out_path, fragment, atomic_symbols, use_published_polarizabilities,use_cm5)

    return charges, pols

def get_chg_pol_from_qchem_output(qchem_out_path, fragment, atomic_symbols, use_published_polarizabilities, use_cm5):
    # parse the output file
    with open(qchem_out_path, "r") as qchem_out:

        # read lines until we read the line before where the volumes are specified
        while True:
            try:
                line = qchem_out.readline()
                if "SCF failed to converge" in line:
                    raise LibraryCallError("Qchem", "charge/polarizability calculation", "SCF failed to converge".format(qchem_out_path), log_path = qchem_out_path)
                line.index("Atom  vol   volFree")
                break
            except ValueError:
                pass

        # read the effective and free volumes from the qchem output into a dictionary that maps from the letter to a list of equivelent values (A -> [0.934, 0.935, 0.934])
        effective_polarizability_dictionary = {}


        # used to keep track of total atoms counted in order to look up atoms in atomic_symbols
        atom_count = 0

        # loop thru all the fragments in the molecule
        at_types = get_atom_types(fragment)

        # loop thru each atom type in the fragment (ie A,B,C in A1B3C2)
        for i in range(0,len(at_types),2):
            atom_index = int(i/2)
            atom_type = at_types[i]

            # if an entry in the effective polarizabilities dictionary does not exist for this atom type, then add a list.
            if atom_type not in effective_polarizability_dictionary:
                effective_polarizability_dictionary[atom_type] = []

            # loop thru each atom of the corresponding atom type (ie 3 times for A3)
            for atom in range(0, int(fragment[atom_index * 2 + 1])):
                # parse the volumes from the next line of the qchem output file
                effective_volume, free_volume = (float(volume) for volume in qchem_out.readline().split()[1:3])
                # look up the free polarizability in the dictionary define at top of file

                if use_published_polarizabilities:
                    free_polarizability = constants.symbol_to_free_polarizability(atomic_symbols[atom_count])
                else:
                    free_polarizability = constants.symbol_to_ccsdt_free_polarizability(atomic_symbols[atom_count])

                # calculate the effective polarizability
                # see equation 6 in https://doi.org/10.1063/1.5009502
                # and equations 5.20 (pg 73) and 5.48 (pg 83) in https://th.fhi-berlin.mpg.de/site/uploads/Publications/Gobre_VV_thesis.pdf 
                effective_polarizability = free_polarizability * (effective_volume / free_volume) ** (4.0/3.0)

                # add the effective polarizability to the corresponding list in the dictionary
                effective_polarizability_dictionary[atom_type].append(effective_polarizability)

                atom_count += 1

        # now construct a list of lists of effective polarizabilities of each atom sorted by fragment, by averaging all the equivelent polarizabilites
        effective_polarizabilities = []

        frag_effective_polarizabilities = []

        # loop thru each atom type in the fragment (ie, A, B, C in A1B3C2)
        for atom_index, atom_type in enumerate(fragment[::2]):

            # loop thru each atom of the corresponding atom symbol (ie 3 times for A3)
            for atom in range(0, int(fragment[atom_index * 2 + 1])):
                # this atom's effective polarizability is the average of all effective polarizabilites for equivelent molecules
                frag_effective_polarizabilities.append(sum(effective_polarizability_dictionary[atom_type]) / len(
                    effective_polarizability_dictionary[atom_type]))

        effective_polarizabilities.append(frag_effective_polarizabilities)

        # read lines until we read the line before where charges are specified
        while True:
            try:
                line = qchem_out.readline()
                if "SCF failed to converge" in line:
                    raise LibraryCallError("Qchem", "charge/polarizability calculation", "SCF failed to converge".format(qchem_out_path), log_path = qchem_out_path)
                if use_cm5:
                  line.index("Charge Model 5")
                else:
                  line.index("Ground-State ChElPG Net Atomic Charges")

                # skip 3 more lines
                for i in range(3):
                    qchem_out.readline()
                break
            except ValueError:
                pass

        # read the charges from the qchem output and construct a dictionary of charges sorted by atom type
        charges_dictionary = {}

        # loop thru each atom type in the fragment (ie A,B,C in A1B3C2)
        for atom_index, atom_type in enumerate(fragment[::2]):

            # if this atom type does not have an entry in the charges dictionary, add one
            if atom_type not in charges_dictionary:
                charges_dictionary[atom_type] = []

            # loop thru each atom of the corresponding atom type (ie 3 times for A3)
            for atom in range(0, int(fragment[atom_index * 2 + 1])):
                # parse the charge from the next line of the qchem output
                charge = float(qchem_out.readline().split()[2])

                charges_dictionary[atom_type].append(charge)

        # now construct list of lists of atom's charges sorted by fragment by averaging all equivelent charges

        charges = []

        frag_charges = []

        # loop thru each atom type in the fragment (ie A,B,C in A1B3C2)
        for atom_index, atom_type in enumerate(fragment[::2]):

            # loop thru each atom of the corresponding atom type (ie 3 times for A3)
            for atom in range(0, int(fragment[atom_index * 2 + 1])):
                # calcualte this atom's charge by averaging all equivelent charges
                frag_charges.append(sum(charges_dictionary[atom_type]) / len(charges_dictionary[atom_type]))

        charges.append(frag_charges)


    return charges, effective_polarizabilities

def get_system_properties(settings_file, config_path, geo_paths,
                          distance_between,
                          use_published_polarizabilities,
                          method,
                          basis,
                          num_digits,
                          virtual_sites,use_cm5):

    # read Settings for the full system
    settings = SettingsReader(settings_file)

    # Obtain the system fragment names
    names = settings.get("molecule", "names")

    print("Generating fitting config file for molecule with fragments: {}.".format(names))

    # Obtain the system properties
    monomer_settings = []
    names = names.split(",")
    fragments = settings.get("molecule", "fragments").split(",")
    charges = settings.get("molecule", "charges").split(",")
    spins = settings.get("molecule", "spins").split(",")
    symmetries = settings.get("molecule", "symmetry").split(",")
    SMILES = settings.get("molecule", "SMILES").split(",")

    # Set up a set of settings for each fragment
    for name, fragment, charge, spin, symmetry, SMILE in zip(names, fragments, charges, spins, symmetries, SMILES):
        monomer_setting = SettingsReader(settings_file)
        monomer_setting.set("molecule", "names", name)
        monomer_setting.set("molecule", "fragments", fragment)
        monomer_setting.set("molecule", "charges", charge)
        monomer_setting.set("molecule", "spins", spin)
        monomer_setting.set("molecule", "symmetry", symmetry)
        monomer_setting.set("molecule", "SMILES", SMILE)
        monomer_settings.append(monomer_setting)

    # Some Sanity checks
    # Check that we have as many geometries as fragments in our molecule
    if len(geo_paths) != len(symmetries):
        raise InconsistentValueError("number of geometries", "number of fragments", len(geo_paths), len(symmetries),
                                     "number of geometries must be equal to the number of fragments in the A3B2_A3B2 type input")

    # split the molecule input string into fragments
    parser = MoleculeInParser("_".join(symmetries))

    fragments = ["".join([atom_type.get_atom_in() for atom_type in frag.get_atom_types()]) for frag in
                 parser.get_fragments()]

    molecule_in = "_".join(
        ["".join([atom_type.get_atom_in() for atom_type in frag.get_atom_and_virtual_site_types()]) for frag in
         parser.get_fragments()])

    chg_list = []
    unrounded_charges = []
    pol_list = []
    c6_list = []

    for i in range(len(symmetries)):
        chgs, pols = calculate_chg_pol_for_config(monomer_settings[i],geo_paths[i], fragments[i], distance_between, use_published_polarizabilities,
                                                     method,
                                                     basis,use_cm5)

        dimer_setting = SettingsReader(settings_file)
        dimer_setting.set("molecule", "names", name + "," + name)
        dimer_setting.set("molecule", "fragments", fragment + "," + fragment)
        dimer_setting.set("molecule", "charges", charge + "," + charge)
        dimer_setting.set("molecule", "spins", spin + "," + spin)
        dimer_setting.set("molecule", "symmetry", symmetry + "," + symmetry)
        dimer_setting.set("molecule", "SMILES", SMILE + "," + SMILE)

        c6 = calculate_c6_for_config(dimer_setting, [monomer_settings[i], monomer_settings[i]], [geo_paths[i],geo_paths[i]],[fragments[i],fragments[i]], distance_between, use_published_polarizabilities,
                                     method,
                                     basis,use_cm5)


        chg_list.append([round(charge, num_digits) for charge in chgs[0]])
        unrounded_charges.append(chgs[0])
        pol_list.append([round(pol, num_digits) for pol in pols[0]])
        c6_list.append([round(c, num_digits) for c in c6[-1]])

    # If we are doing a dimer, we need the intermolecular c6
    if len(symmetries) == 2:
        true_dimer_setting = SettingsReader(settings_file)
        c6 = calculate_c6_for_config(true_dimer_setting, [monomer_settings[0], monomer_settings[1]], [geo_paths[0],geo_paths[1]],[fragments[0],fragments[1]], distance_between, use_published_polarizabilities,
                                     method,
                                     basis,use_cm5)

        c6_list.append([round(c, num_digits) for c in c6[-1]])

    system.format_print("Unrounded charges: {}".format(unrounded_charges),
                        italics=True)
    system.format_print("Rounded charges: {}".format(chg_list),
                        italics=True)

    expected_charge = sum([int(charge) for charge in charges])
    unrounded_charge = sum([charge for sub_list in unrounded_charges for charge in sub_list])
    rounded_charge = sum([charge for sub_list in chg_list for charge in sub_list])

    system.format_print("Expected charge: {}; Unrounded charge: {}; Rounded charge: {}".format(expected_charge,

 unrounded_charge,

 rounded_charge),
                        italics=True)

    system.format_print("Difference between expected and unrounded charge: {}.".format(expected_charge - unrounded_charge),
                        italics=True)
    system.format_print("Difference between expected and rounded charge: {}.".format(expected_charge - rounded_charge),
                        italics=True)

    return chg_list, pol_list, c6_list[-1]

def write_config_file(settings_file, config_path, chg, 
                      pols, geo_paths, C6, polfacs, 
                      d6, A, 
                      kmin, kmax, dmin, dmax, bmin, bmax,
                      kmin_init, kmax_init, dmin_init, dmax_init,
                      bmin_init, bmax_init, r_in, r_out,
                      energy_range, alpha,
                      virtual_sites_label,
                      var_intra, var_inter, var_virtual_sites):

    settings = SettingsReader(settings_file)

    # Obtain the system properties
    monomer_settings = []
    names = settings.get("molecule","names").split(",")
    fragments = settings.get("molecule", "fragments").split(",")
    charges = settings.get("molecule", "charges").split(",")
    spins = settings.get("molecule", "spins").split(",")
    symmetries = settings.get("molecule", "symmetry").split(",")
    SMILES = settings.get("molecule", "SMILES").split(",")

    # Set up a set of settings for each fragment
    for name, fragment, charge, spin, symmetry, SMILE in zip(names, fragments, charges, spins, symmetries, SMILES):
        monomer_setting = SettingsReader(settings_file)
        monomer_setting.set("molecule", "names", name)
        monomer_setting.set("molecule", "fragments", fragment)
        monomer_setting.set("molecule", "charges", charge)
        monomer_setting.set("molecule", "spins", spin)
        monomer_setting.set("molecule", "symmetry", symmetry)
        monomer_setting.set("molecule", "SMILES", SMILE)
        monomer_settings.append(monomer_setting)


    molecule_in = "_".join(symmetries)
    parser = MoleculeInParser(molecule_in)

    fragments = []

    for geo_path, setting in zip(geo_paths, monomer_settings):
        with open(geo_path, "r") as geo_file:
            frag_string = "\n".join(geo_file.read().splitlines()[2:])
            fragments.append(
                Fragment.read_xyz(frag_string, setting.get("molecule", "names"), setting.getint("molecule", "charges"),
                                  setting.getint("molecule", "spins"), setting.get("molecule", "SMILES"),
                                  setting.get("molecule", "symmetry")))    

    molecule = Molecule(fragments)

    # create the config file!
    configwriter = configparser.ConfigParser()
    configwriter.add_section("common")
    configwriter.add_section("fitting")

    configwriter.set("common", "molecule", molecule_in)

    configwriter.set("fitting", "number_of_atoms", str(molecule.get_num_atoms()))
    configwriter.set("fitting", "number_of_electrostatic_sites", str(molecule.get_num_atoms() + sum([int(i) for i in settings.get("molecule", "use_mbpol").split(",")])))


    excluded_pairs12, excluded_pairs13, excluded_pairs14 = molecule.get_excluded_pairs()

    configwriter.set("fitting", "excluded_pairs_12", "{}".format(excluded_pairs12))
    configwriter.set("fitting", "excluded_pairs_13", "{}".format(excluded_pairs13))
    configwriter.set("fitting", "excluded_pairs_14", "{}".format(excluded_pairs14))

    configwriter.set("fitting", "charges", str(chg))
    configwriter.set("fitting", "polarizabilities", str(pols))
    if polfacs is None:
        configwriter.set("fitting", "polarizability_factors", str(pols))
    else:
        configwriter.set("fitting", "polarizability_factors", str(polfacs))

    configwriter.set("fitting", "k_min", str(kmin))
    configwriter.set("fitting", "k_max", str(kmax))
    configwriter.set("fitting", "d_min", str(dmin))
    configwriter.set("fitting", "d_max", str(dmax))

    configwriter.set("fitting", "k_min_init", str(kmin_init))
    configwriter.set("fitting", "k_max_init", str(kmax_init))
    configwriter.set("fitting", "d_min_init", str(dmin_init))
    configwriter.set("fitting", "d_max_init", str(dmax_init))

    configwriter.set("fitting", "b_min", str(bmin))
    configwriter.set("fitting", "b_max", str(bmax))
    configwriter.set("fitting", "b_min_init", str(bmin_init))
    configwriter.set("fitting", "b_max_init", str(bmax_init))

    configwriter.set("fitting", "r_in", str(r_in))
    configwriter.set("fitting", "r_out", str(r_out))

    configwriter.set("fitting", "C6", str(C6))
    if d6 is None:
        configwriter.set("fitting", "d6", str([0.0]*len(C6)))
    else:
        configwriter.set("fitting", "d6", str(d6))
    if A is None:
        configwriter.set("fitting", "A", str([0.0]*len(C6)))
    else:
        configwriter.set("fitting", "A", str(A))

    configwriter.set("fitting", "var_intra", var_intra)
    configwriter.set("fitting", "var_inter", var_inter)
    configwriter.set("fitting", "var_virtual_sites", var_virtual_sites)
    configwriter.set("fitting", "alpha", str(alpha))
    configwriter.set("fitting", "energy_range", str(energy_range))
    configwriter.set("fitting", "virtual_site_labels", "{}".format(virtual_sites_label).replace("'","").replace('"',""))

    config_path = files.init_file(config_path, files.OverwriteMethod.get_from_settings(settings))

    with open(config_path, "w") as config_file:
        configwriter.write(config_file)

    print("Completed generating config file {}.".format(config_path))

