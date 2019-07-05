import os
import math, configparser
from collections import OrderedDict

from potential_fitting.utils import constants, SettingsReader, files
from potential_fitting.exceptions import InvalidValueError, InconsistentValueError
from potential_fitting.molecule import Molecule, Fragment, xyz_to_molecules
from potential_fitting.polynomials import MoleculeInParser

qchem_template = "qchem_template"

def make_config(settings_file, molecule_in, config_path, *geo_paths, distance_between = 20, use_published_polarizabilities = True):
    """
    Generates the config file for the fitcode for the given geometries

    Args:
        settings_file - the file containing relevent settings information
        molecule_in - A3B1 formatted string
        config_path - path to file to write config file to, should end in .ini
        geo_paths - paths to each geometry to include in the config, should be 1 to 3 of them (inclusive)
        distance_between - the distance between each geometry, in angstroms
        use_published_polarizabilities - use the polarizabilities from the 2018 Schwerdtfeger & Nagle paper; otherwise, use those calculated using MolPRO with ccsd(t)
    """


    settings = SettingsReader(settings_file)

    print("Generating fitting config file for molecule with fragments: {}".format(settings.get("molecule", "names")))

    print("Preparing qchem input...")

    monomer_settings = []
    names = settings.get("molecule", "names").split(",")
    fragments = settings.get("molecule", "fragments").split(",")
    charges = settings.get("molecule", "charges").split(",")
    spins = settings.get("molecule", "spins").split(",")
    symmetries = settings.get("molecule", "symmetry").split(",")


    for name, fragment, charge, spin, symmetry in zip(names, fragments, charges, spins, symmetries):
        monomer_setting = SettingsReader(settings_file)
        monomer_setting.set("molecule", "names", name)
        monomer_setting.set("molecule", "fragments", fragment)
        monomer_setting.set("molecule", "charges", charge)
        monomer_setting.set("molecule", "spins", spin)
        monomer_setting.set("molecule", "symmetry", symmetry)
        monomer_settings.append(monomer_setting)

    # split the molecule input string into fragments

    parser = MoleculeInParser(molecule_in)

    fragments = ["".join([atom_type.get_atom_in() for atom_type in frag.get_atom_types()]) for frag in parser.get_fragments()]

    molecule_in = "_".join(fragments)

    if len(geo_paths) != len(fragments):
        raise InconsistentValueError("number of geometries", "number of fragments", len(geo_paths), len(fragments), "number of geometries must be equal to the number of fragments in the A3B2_A3B2 type input")

    qchem_in_path = os.path.join(settings.get("files", "log_path"), "get_config_qchem.in")
    with open(qchem_in_path, "w") as qchem_in:
        
        # tells qchem that the molecule has started
        qchem_in.write("$molecule\n")

        # charge and spin line
        qchem_in.write("{} {}\n".format(sum([int(charge) for charge in settings.get("molecule", "charges").split(",")]), 1 + sum([int(spin) - 1 for spin in settings.get("molecule", "spins").split(",")])))

        atomic_symbols = []

        # if there are 0 geometries specified, raise an error
        if len(geo_paths) == 0:
            raise InvalidValueError("number of geometries", len(geo_paths), "at least 1")

        # if there is at least 1 geometry specified
        else:

            # read geometry into a molecule object
            molecule1 = xyz_to_molecules(geo_paths[0], settings = monomer_settings[0])[0]

            # move molecule1 to its standard orientation
            molecule1.move_to_center_of_mass()
            molecule1.rotate_on_principal_axes()

            # copy this molecule to the qchem input
            qchem_in.write(molecule1.to_xyz())
            qchem_in.write("\n")

            # add molecule1's atoms to the list of atomic_symbols
            for atom in molecule1.get_atoms():
                atomic_symbols.append(atom.get_name())

            # if there are at least 2 geometries specified
            if len(geo_paths) > 1:

                # read geonetry into molecule object
                molecule2 = xyz_to_molecules(geo_paths[1], settings = monomer_settings[1])[0]

                # move molecule2 to its standard orientation
                molecule2.move_to_center_of_mass()
                molecule2.rotate_on_principal_axes()

                # move this molecule away from the first one by distance_between angstroms
                molecule2.translate(distance_between, 0, 0)

                # copy this molecule to the qchem input
                qchem_in.write(molecule2.to_xyz())
                qchem_in.write("\n")

                # add molecule2's atoms to the list of atomic_symbols
                for atom in molecule2.get_atoms():
                    atomic_symbols.append(atom.get_name())

                # if there are at least 3 geometries specified
                if len(geo_paths) > 2:

                    # read geometry into a molecule object
                    molecule3 = xyz_to_molecules(geo_paths[2], settings = monomer_settings[2])[0]

                    # move molecule3 to its standard orientation
                    molecule3.move_to_center_of_mass()
                    molecule3.rotate_on_principal_axes()

                    # move molecule3 so it is equadist from the other two molecules
                    molecule3.translate(distance_between/2, distance_between * math.sqrt(3) / 2, 0)

                    # copy this molecule to the qchem input
                    qchem_in.write(molecule3.to_xyz())
                    qchem_in.write("\n")

                    # add molecule3's atoms to the list of atomic_symbols
                    for atom in molecule3.get_atoms():
                        atomic_symbols.append(atom.get_name())
                

        # tells qchem that the molecule has ended
        qchem_in.write("$end\n")

        qchem_in.write("$rem\n")
        qchem_in.write("method " + settings.get("config", "method", "wb97m-v") + "\n")
        qchem_in.write("basis " + settings.get("config", "basis", "aug-cc-pvtz") + "\n")

        # read the qchem template and append it to the qchem in
        qchem_template_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), qchem_template)
        with open(qchem_template_path, "r") as template:
            for line in template:
                qchem_in.write(line)

    qchem_out_path = os.path.join(settings.get("files", "log_path"), "get_config_qchem.out")
    qchem_log_path = os.path.join(settings.get("files", "log_path"), "get_config_qchem.log")

    num_threads = settings.getint("qchem", "num_threads")

    print("Executing qchem calculation...")
    # perform qchem system call
    #os.system("qchem -nt {} {} {} > {}".format(num_threads, qchem_in_path, qchem_out_path, qchem_log_path))

    print("Parsing qchem output...")

    # parse the output file
    with open(qchem_out_path, "r") as qchem_out:
        
        # read lines until we read the line before where the volumes are specified
        while True:
            try:
                line = qchem_out.readline()
                if line == "":
                    print("Something went wrong with qchem")
                    exit(1)
                line.index("Atom  vol   volFree")
                break
            except ValueError:
                pass

        # read the effective and free volumes from the qchem output into a dictionary that maps from the letter to a list of equivelent values (A -> [0.934, 0.935, 0.934])
        effective_polarizability_dictionary = {}

        # used to keep track of total atoms counted in order to look up atoms in atomic_symbols
        atom_count = 0

        # loop thru all the fragments in the molecule
        for fragment in fragments:

            # loop thru each atom type in the fragment (ie A,B,C in A1B3C2)
            for atom_index, atom_type in enumerate(fragment[::2]):

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
                    effective_polarizability = free_polarizability * effective_volume / free_volume

                    # add the effective polarizability to the corresponding list in the dictionary
                    effective_polarizability_dictionary[atom_type].append(effective_polarizability)

                    atom_count += 1

        # now construct a list of lists of effective polarizabilities of each atom sorted by fragment, by averaging all the equivelent polarizabilites
        effective_polarizabilities = []

        # loop over all the fragments in the molecule
        for fragment in fragments:

            frag_effective_polarizabilities = []

            # loop thru each atom type in the fragment (ie, A, B, C in A1B3C2)
            for atom_index, atom_type in enumerate(fragment[::2]):

                # loop thru each atom of the corresponding atom symbol (ie 3 times for A3)
                for atom in range(0, int(fragment[atom_index * 2 + 1])):

                    # this atom's effective polarizability is the average of all effective polarizabilites for equivelent molecules
                    frag_effective_polarizabilities.append(sum(effective_polarizability_dictionary[atom_type]) / len(effective_polarizability_dictionary[atom_type]))

            effective_polarizabilities.append(frag_effective_polarizabilities)
        
        # read lines until we read the line before where c6 constants are specified
        while True:
            try:
                qchem_out.readline().index("D2X   (i,j)")
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

        while(True):

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

            while(True):

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

        # read lines until we read the line before where charges are specified
        while True:
            try:
                qchem_out.readline().index("Ground-State ChElPG Net Atomic Charges")
                # skip 3 more lines
                for i in range(3):
                    qchem_out.readline()
                break
            except ValueError:
                pass

        # read the charges from the qchem output and construct a dictionary of charges sorted by atom type
        charges_dictionary = {}

        # loop thru all the fragments in the molecule
        for fragment in fragments:

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

        # loop thru all the fragments in the molecule
        for fragment in fragments:

            frag_charges = []

            # loop thru each atom type in the fragment (ie A,B,C in A1B3C2)
            for atom_index, atom_type in enumerate(fragment[::2]):

                # loop thru each atom of the corresponding atom type (ie 3 times for A3)
                for atom in range(0, int(fragment[atom_index * 2 + 1])):

                    # calcualte this atom's charge by averaging all equivelent charges
                    frag_charges.append(sum(charges_dictionary[atom_type]) / len(charges_dictionary[atom_type]))

            charges.append(frag_charges)

    print("Writing config file...")

    # create the config file!
    configwriter = configparser.ConfigParser()
    configwriter.add_section("common")
    configwriter.add_section("fitting")

    configwriter.set("common", "molecule", molecule_in)

    configwriter.set("fitting", "number_of_atoms", str(len(atomic_symbols)))
    configwriter.set("fitting", "number_of_electrostatic_sites", str(len(atomic_symbols)))

    fragments = []

    for geo_path, setting in zip(geo_paths, monomer_settings):
        with open(geo_path, "r") as geo_file:
            frag_string = "\n".join(geo_file.read().splitlines()[2:])
            fragments.append(Fragment.read_xyz(frag_string, setting.get("molecule", "names"), setting.getint("molecule", "charges"), setting.getint("molecule", "spins"), setting.get("molecule", "SMILES"), setting.get("molecule", "symmetry")))

    molecule = Molecule(fragments)

    excluded_pairs12, excluded_pairs13, excluded_pairs14 = molecule.get_excluded_pairs()

    configwriter.set("fitting", "excluded_pairs_12", "{}".format(excluded_pairs12))
    configwriter.set("fitting", "excluded_pairs_13", "{}".format(excluded_pairs13))
    configwriter.set("fitting", "excluded_pairs_14", "{}".format(excluded_pairs14))

    configwriter.set("fitting", "charges", str(charges))
    configwriter.set("fitting", "polarizabilities", str(effective_polarizabilities))
    configwriter.set("fitting", "polarizability_fractions", str(effective_polarizabilities))

    configwriter.set("fitting", "k_min", str(0.0))
    configwriter.set("fitting", "k_max", str(5.0))
    configwriter.set("fitting", "d_min", str(0.0))
    configwriter.set("fitting", "d_max", str(7.0))

    configwriter.set("fitting", "C6", str(c6_constants))
    configwriter.set("fitting", "d6", str([[0 for x in c6] for c6 in c6_constants]))

    configwriter.set("fitting", "var_intra", "exp")
    configwriter.set("fitting", "var_inter", "exp")
    configwriter.set("fitting", "var_lonepairs", "coul")
    configwriter.set("fitting", "energy_range", str(50.0))
    configwriter.set("fitting", "virtual_site_labels", "[X,Y,Z]")

    config_path = files.init_file(config_path, files.OverwriteMethod.get_from_settings(settings))

    with open(config_path, "w") as config_file:
        configwriter.write(config_file)

    print("Completed generating config file {}.".format(config_path))
