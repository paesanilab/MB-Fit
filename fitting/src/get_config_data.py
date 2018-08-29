import sys, os
import math, configparser
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)) + "/../../")
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)) + "/../../qm_mb_energy_calculator/src")

import constants, settings_reader
from collections import OrderedDict
from exceptions import InvalidValueError, InconsistentValueError
import molecule_parser

qchem_template = "/qchem_template"

# table of free polarizabilities for each atom
free_polarizabilities = {
    "H":    0.66582,
    "B":    3.64084,
    "C":    1.84613,
    "N":    1.08053,
    "O":    0.86381,
    "F":    0.50682,
    "P":    3.72507,
    "S":    3.24768
}

def make_config(settings_file, molecule_in, config_path, *geo_paths, distance_between = 20):
    """
    Generates the config file for the fitcode for the given geometries

    Args:
        settings_file - the file containint relevent settings information
        molecule_in - A3B1 formatted string
        config_path - path to file to write config file to, should end in .ini
        geo_paths - paths to each geometry to include in the config, should be 1 to 3 of them (inclusive)
        distance_between - the distance between each geometry, in angstroms
    """

    settings = settings_reader.SettingsReader(settings_file)

    # split the molecule input string into fragments
    fragments = molecule_in.split("_")

    if len(geo_paths) != len(fragments):
        raise InconsistentValueError("number of geometries", "number of fragments", len(geo_paths), len(fragments), "number of geometries must be equal to the number of fragments in the A3B2_A3B2 type input")

    log_directory = settings.get("files", "log_path")

    with open(log_directory + "/qchem.in", "w") as qchem_in:
        
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
            molecule1 = molecule_parser.xyz_to_molecules(geo_paths[0])[0]

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
                molecule2 = molecule_parser.xyz_to_molecules(geo_paths[1])[0]

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
                    molecule3 = molecule_parser.xyz_to_molecules(geo_paths[2])[0]

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

        # read the qchem template and append it to the qchem in
        with open(os.path.dirname(os.path.abspath(__file__)) + "/" + qchem_template, "r") as template:
            for line in template:
                qchem_in.write(line)

    # perform qchem system call
    #os.system("qchem -nt 4 {} {} > {}".format(log_directory + "/qchem.in", log_directory + "/qchem.out", log_directory + "/qchem.log"))

    # parse the output file
    with open(log_directory + "/qchem.out", "r") as qchem_out:
        
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

            # loop thru each atomic symbol in the fragment (ie A,B,C in A1B3C2)
            for atom_index, atom_type in enumerate(fragment[::2]):

                # if an entry in the effective polarizabilities dictionary does not exist for this atom type, then add a list.
                if atom_type not in effective_polarizability_dictionary:
                    effective_polarizability_dictionary[atom_type] = []

                # loop thru each atom of the corresponding atomic symbol (ie 3 times for A3)
                for atom in range(0, int(fragment[atom_index * 2 + 1])):

                    # parse the volumes from the next line of the qchem output file
                    effective_volume, free_volume = (float(volume) for volume in qchem_out.readline().split()[1:3])
                    # look up the free polarizability in the dictionary define at top of file
                    free_polarizability = free_polarizabilities[atomic_symbols[atom_count]]

                    # calculate the effective polarizability
                    effective_polarizability = free_polarizability * effective_volume / free_volume

                    # add the effective polarizability to the corresponding list in the dictionary
                    effective_polarizability_dictionary[atom_type].append(effective_polarizability)

                    atom_count += 1

        # now construct a list of lists of effective polarizabilities of each atom sorted by fragment
        effective_polarizabilities = []

        # loop over all the fragments in the molecule
        for fragment in fragments:

            frag_effective_polarizabilities = []

            # loop thru each atomic symbol in the fragment (ie, A, B, C in A1B3C2)
            for atom_index, atom_type in enumerate(fragment[::2]):

                # loop thru each atom of the corresponding atomic symbol (ie 3 times for A3)
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

        c6_constant_lists = [OrderedDict() for x in range(len(fragments) + 1)]

        # read the c6 constants of each pair of atoms

        fragment_index_a = 0
        atom_index_a = 0
        atom_a = 0

        while(True):

            if atom_a >= int(fragments[fragment_index_a][atom_index_a * 2 + 1]):
                atom_index_a += 1
                atom_a = 0

                if atom_index_a >= len(fragments[fragment_index_a][::2]):
                    fragment_index_a += 1
                    atom_index_a = 0

                    if fragment_index_a >= len(fragments):
                        break

            fragment_index_b = fragment_index_a
            atom_index_b = atom_index_a
            atom_b = atom_a + 1

            while(True):
                if atom_b >= int(fragments[fragment_index_b][atom_index_b * 2 + 1]):
                    atom_index_b += 1
                    atom_b = 0

                    if atom_index_b >= len(fragments[fragment_index_b][::2]):
                        fragment_index_b += 1
                        atom_index_b = 0

                        if fragment_index_b >= len(fragments):
                            break

                line = qchem_out.readline()

                pair = [int(index) for index in line.split()[1].replace("(", "").replace(")", "").split(",")]

                c6 = float(line.split()[2])

                # convert c6 constant from Haretrees/Bohr^6 to kcal/mol/A^6
                c6 *= constants.au_per_bohr6_to_kcal_per_ang6

                atomic_symbol_a = fragments[fragment_index_a][::2][atom_index_a]
                atomic_symbol_b = fragments[fragment_index_b][::2][atom_index_b]
                
                # check if this is an intrafragmental c6 constant
                if fragment_index_a == fragment_index_b:
                    if atomic_symbol_a + atomic_symbol_b in c6_constant_lists[fragment_index_a]:
                        c6_constant_lists[fragment_index_a][atomic_symbol_a + atomic_symbol_b].append(c6)

                    elif atomic_symbol_b + atomic_symbol_a in c6_constant_lists[fragment_index_a]:
                        c6_constant_lists[fragment_index_a][atomic_symbol_b + atomic_symbol_a].append(c6)

                    else:
                        c6_constant_lists[fragment_index_a][atomic_symbol_a + atomic_symbol_b] = []
                        c6_constant_lists[fragment_index_a][atomic_symbol_a + atomic_symbol_b].append(c6)
                # otherwise this is an interfragmental c6 constant
                else:
                    if atomic_symbol_a + atomic_symbol_b in c6_constant_lists[len(fragments)]:
                        c6_constant_lists[len(fragments)][atomic_symbol_a + atomic_symbol_b].append(c6)

                    elif atomic_symbol_b + atomic_symbol_a in c6_constant_lists[len(fragments)]:
                        c6_constant_lists[len(fragments)][atomic_symbol_b + atomic_symbol_a].append(c6)

                    else:
                        c6_constant_lists[len(fragments)][atomic_symbol_a + atomic_symbol_b] = []
                        c6_constant_lists[len(fragments)][atomic_symbol_a + atomic_symbol_b].append(c6)

                atom_b += 1

            atom_a += 1

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

            frag_charges = []

            # loop thru each atomic symbol in the fragment (ie A,B,C in A1B3C2)
            for atom_index, atom_type in enumerate(fragment[::2]):

                # if this atom type does not have an entry in the charges dictionary, add one
                if atom_type not in charges_dictionary:
                    charges_dictionary[atom_type] = []

                # loop thru each atom of the corresponding atomic symbol (ie 3 times for A3)
                for atom in range(0, int(fragment[atom_index * 2 + 1])):

                    # parse the charge from the next line of the qchem output
                    charge = float(qchem_out.readline().split()[2])

                    charges_dictionary[atom_type].append(charge)

        # now construct list of lists of atom's charges sorted by fragment by averaging all equivelent charges

        charges = []

        # loop thru all the fragments in the molecule
        for fragment in fragments:

            frag_charges = []

            # loop thru each atomic symbol in the fragment (ie A,B,C in A1B3C2)
            for atom_index, atom_type in enumerate(fragment[::2]):

                # loop thru each atom of the corresponding atomic symbol (ie 3 times for A3)
                for atom in range(0, int(fragment[atom_index * 2 + 1])):

                    # calcualte this atom's charge by averaging all equivelent charges
                    frag_charges.append(sum(charges_dictionary[atom_type]) / len(charges_dictionary[atom_type]))

            charges.append(frag_charges)

    # create the config file!
    configwriter = configparser.ConfigParser()
    configwriter.add_section("common")
    configwriter.add_section("fitting")

    configwriter.set("common", "molecule", molecule_in)
    configwriter.set("common", "polynomial_order", settings.get("poly_generator", "order"))

    configwriter.set("fitting", "number_of_atoms", str(len(atomic_symbols)))
    configwriter.set("fitting", "number_of_electrostatic_sites", str(len(atomic_symbols)))

    configwriter.set("fitting", "excluded_pairs_12", "[TODO: User must fill this in]")
    configwriter.set("fitting", "excluded_pairs_13", "[TODO: User must fill this in]")
    configwriter.set("fitting", "excluded_pairs_14", "[TODO: User must fill this in]")

    configwriter.set("fitting", "charges", str(charges))
    configwriter.set("fitting", "polarizabilities", str(effective_polarizabilities))
    configwriter.set("fitting", "polarizability_fractions", str(effective_polarizabilities))

    configwriter.set("fitting", "k_min", str(0.0))
    configwriter.set("fitting", "k_max", str(5.0))
    configwriter.set("fitting", "d_min", str(0.0))
    configwriter.set("fitting", "d_max", str(7.0))

    configwriter.set("fitting", "C6", str(c6_constants))
    configwriter.set("fitting", "d6", str([[0 for x in c6] for c6 in c6_constants]))

    configwriter.set("fitting", "var", "exp")
    configwriter.set("fitting", "energy_range", str(50.0))
    configwriter.set("fitting", "virtual_site_labels", str(["X", "Y", "Z"]))

    with open(config_path, "w") as config_file:
        configwriter.write(config_file)
