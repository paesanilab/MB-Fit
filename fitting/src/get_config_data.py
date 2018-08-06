import sys
sys.path.insert(0, "../../")

import os
import constants
import configparser

qchem_template = "./qchem_template"

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

def get_config_data(settings, geo_path):

    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)

    # read number of atoms per fragment from config
    atoms_per_frag = [int(num_atoms) for num_atoms in config["molecule"]["fragments"].split(",")]

    # find the total number of ragments
    number_of_fragments = len(atoms_per_frag)

    # find the total number of atoms by summing the number of atoms in each fragment
    number_of_atoms = sum(atoms_per_frag)

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
        effective_volumes = []
        free_volumes = []
        for i in range(number_of_atoms):
            effective_volume, free_volume = (float(volume) for volume in qchem_out.readline().split()[1:3])
            effective_volumes.append(effective_volume)
            free_volumes.append(free_volume)

        # look up the free polarizabilities of each atom
        free_polars = [free_polarizabilities[atomic_symbol] for atomic_symbol in atomic_symbols]

        # calculate the effective polarizability of each atom
        effective_polars = [(free_polar * effective_volume / free_volume) for free_polar, effective_volume, free_volume in zip(free_polars, effective_volumes, free_volumes)]

        # TODO: AVERAGE EQUIVELENT POLARIZABILITIES
        
        # read lines until we read the line before where c6 constants are specified
        while True:
            try:
                qchem_out.readline().index("D2X   (i,j)")
                break
            except ValueError:
                pass

        c6_constants = [[] for x in range(number_of_fragments + 1)]
        
        # read the c6 constants of each pair of atoms
        for i in range(sum(range(number_of_atoms))): # this is how many lines will be devoted to c6 constants in qchem output
            line = qchem_out.readline()
            pair = [int(index) for index in line.split()[1].replace("(", "").replace(")", "").split(",")]

            c6 = float(line.split()[2])

            # convert c6 constant from Haretrees/Bohr^6 to kcal/mol/A^6
            c6 *= constants.au_per_bohr6_to_kcal_per_ang6
            
            # check if this pair is intra-molecular
            intra_molecular = False
            for i in range(number_of_fragments):
                if pair[0] >= sum(atoms_per_frag[:i]) and pair[1] >= sum(atoms_per_frag[:i]) and pair[0] < sum(atoms_per_frag[:i + 1]) and pair[1] < sum(atoms_per_frag[:i + 1]):
                    intra_molecular = True
                    c6_constants[i].append(c6)
                    break
            
            if not intra_molecular:
                c6_constants[number_of_fragments].append(c6)

        # TODO: AVERAGE EQUIVELENT C6 CONSTANTS

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

        charges = []

        # read the charges for each atom
        for i in range(number_of_atoms):
            charge = float(qchem_out.readline().split()[2])
            charges.append(charge)

        # TODO: AVERAGE EQUIVELENT CHARGES

        print("Charges: {}".format(charges))
        print("Polars: {}".format(effective_polars))
        print("C6: {}".format(c6_constants))
        return charges, effective_polars, c6_constants

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python get_config_data.py <settings> <optimized geometry>")
        exit(1)
    get_config_data(sys.argv[1], sys.argv[2])
