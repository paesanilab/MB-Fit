import sys, os
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)) + "/../../")

import constants
import configparser
from collections import OrderedDict

qchem_template = os.path.dirname(__file__) + "/qchem_template"

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

def make_config(settings, molecule_in, geo_path, config_path):

    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)

    # split the molecule input string into fragments
    fragments = molecule_in.split("_")

    log_directory = config["files"]["log_path"]

    with open(log_directory + "/qchem.in", "w") as qchem_in:
        
        # tells qchem that the molecule has started
        qchem_in.write("$molecule\n")

        # charge and spin line
        qchem_in.write("{} {}\n".format(sum([int(charge) for charge in config["molecule"]["charges"].split(",")]), 1 + sum([int(spin) - 1 for spin in config["molecule"]["spins"].split(",")])))

        # tells qchem that the molecule has ended

        # read the geometry from the input xyz file and write it to the qchem input file
        with open(geo_path, "r") as geo_file:
            # skip atom count and comment lines
            geo_file.readline()
            geo_file.readline()
            

            atomic_symbols = []
            # copy lines over from input geometry to qchem in
            for line in geo_file.readlines():
                atomic_symbols.append(line.split()[0])
                qchem_in.write(line)

        qchem_in.write("$end\n")

        # read the qchem template and append it to the qchem in
        with open(qchem_template, "r") as template:
            for line in template:
                qchem_in.write(line)

    # the qchem input file should now be generated

    # perform qchem system call
    os.system("qchem -nt 4 {} {} > {}".format(log_directory + "/qchem.in", log_directory + "/qchem.out", log_directory + "/qchem.log"))

    # parse the output file
    with open(log_directory + "/qchem.out", "r") as qchem_out:
        
        # read lines until we read the line before where the volumes are specified
        while True:
            try:
                qchem_out.readline().index("Atom  vol   volFree")
                break
            except ValueError:
                pass

        # read the effective and free volumes from the qchem output
        effective_polarizabilities = []
        atom_count = 0

        # loop thru all the fragments in the molecule
        for fragment in fragments:

            frag_effective_polarizabilities = []

            # loop thru each atomic symbol in the fragment (ie A,B,C in A1B3C2)
            for atom_index, atom_type in enumerate(fragment[::2]):

                # loop thru each atom of the corresponding atomic symbol (ie 3 times for A3)
                for atom in range(0, int(fragment[atom_index * 2 + 1])):

                    # parse the volumes from the next line of the qchem output file
                    effective_volume, free_volume = (float(volume) for volume in qchem_out.readline().split()[1:3])
                    # look up the free polarizability in the dictionary define at top of file
                    free_polarizability = free_polarizabilities[atomic_symbols[atom_count]]

                    # calculate the effective polarizability
                    effective_polarizability = free_polarizability * effective_volume / free_volume

                    frag_effective_polarizabilities.append(effective_polarizability)

                    atom_count += 1
            
            effective_polarizabilities.append(frag_effective_polarizabilities)

        # TODO: AVERAGE EQUIVELENT POLARIZABILITIES
        
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

        # SAVED FOR LATER!    charge = float(qchem_out.readline().split()[2])

        # read the charges from the qchem output
        charges = []
        atom_count = 0

        # loop thru all the fragments in the molecule
        for fragment in fragments:

            frag_charges = []

            # loop thru each atomic symbol in the fragment (ie A,B,C in A1B3C2)
            for atom_index, atom_type in enumerate(fragment[::2]):

                # loop thru each atom of the corresponding atomic symbol (ie 3 times for A3)
                for atom in range(0, int(fragment[atom_index * 2 + 1])):

                    # parse the charge from the next line of the qchem output
                    charge = float(qchem_out.readline().split()[2])
                    frag_charges.append(charge)

                    atom_count += 1
            
            charges.append(frag_charges)

        # TODO: Average charges

    # create the config file!
    configwriter = configparser.ConfigParser()
    configwriter.add_section("common")
    configwriter.add_section("fitting")

    configwriter.set("common", "molecule", molecule_in)

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

    configwriter.set("fitting", "C6", str(c6_constants[len(fragments)]))
    configwriter.set("fitting", "d6", str([0 for x in c6_constants[len(fragments)]]))

    configwriter.set("fitting", "var", "exp")
    configwriter.set("fitting", "energy_range", str(50.0))
    configwriter.set("fitting", "virtual_site_lables", str(["X", "Y", "Z"]))

    with open(config_path, "w") as config_file:
        configwriter.write(config_file)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python get_config_data.py <settings> <A1B2> <optimized geometry> <config_path>")
        exit(1)
    make_config(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
