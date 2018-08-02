import os
import constants
import configparser

qchem_template = "./qchem_template"

free_polarizabilities = {
                        
}

def get_config_data(settings, geo_path, log_directory):

    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)

    # read number of atoms per fragment from config
    atoms_per_frag = int(num_atoms) for num_atoms in config["molecule"]["fragments"].split(",")

    # find the total number of atoms by summing the number of atoms in each fragment
    number_of_atoms = sum(atoms_per_frag)

    with open(log_directory + "/qchem.in", "w") as qchem_in:
        
        # tells qchem that the molecule has started
        qchem_in.write("$molecule\n")

        # read the geometry from the input xyz file and write it to the qchem input file
        # charge and multiplicity must be in comment line
        with open(geo_path, "r") as geo_file:
            atomic_symbols = []
            # copy lines over from input geometry to qchem in
            for line in geo_file.readlines():
                atomic_symbols.append(line.split()[0])
                qchem_in.write(line)

        # tells qchem that the molecule has ended
        qchem_in.write("$end\n")

        # read the qchem template and append it to the qchem in
        with open(qchem_template, "r") as template:
            for line in template:
                qchem_in.write(line)

    # the qchem input file should now be generated

    # perform qchem system call

    os.system("qchem -nt 4 {} {} > {}".format(log_directory + "/qchem.in", log_directory + "/qchem.out", log_directory + "/qchem.log"))
    

    # parese the output file
    with open(log_directory + "/qchem.out", "r") as qchem_out:
        
        # read lines until we read the line before where the volumes are specified
        while True:
            try:
                qchem_out.readline().index("Atom vol volFree")
                break
            except IndexError:
                pass

        # read the effective and free volumes from the qchem output
        effective_volumes = []
        free_volumes = []
        for i in range(number_of_atoms):
            effective_volume, free_volume = float(volume) for volume in qchem_out.readline().split()[1:3]
            effective_volumes.append(effective_volume)
            free_volumes.append(free_volumes)

        # look up the free polarizabilities of each atom
        free_polars = free_polarizabilities[atomic_symbol] for atomic_symbol in atomic_symbols

        # calculate the effective polarizability of each atom
        effective_polars = free_polar * effective_volume / free_volume for free_polar, effective_volume, free_volume in zip(free_polar, effective_volumes, free_volumes)

        # TODO: AVERAGE EQUIVELENT POLARIZABILITIES
        
        # read lines until we read the line before where c6 constants are specified
        while True:
            try:
                qchem_out.readline().index("D2X   (i,j)")
                break
            except IndexError:
                pass

        c6_constants = []
        
        # read the c6 constants of each pair of atoms
        for i in range(sum(range(number_of_atoms))): # this is how many lines will be devoted to c6 constants in qchem output
            line = qchem_out.readline()
            pair = int(index) for index in line.split()[1].replace("(", "").replace(")", "").split(",")

            # check if this pair is intra-molecular
            intra_molecular = False
            for i in range(atoms_per_frag):
                if pair[0] > sum(atoms_per_frag[:i]) and pair[1] > sum(atoms_per_frag[:i]) and pair[0] < sum(atoms_per_frag[:i + 1]) and pair[1] > sum(atoms_per_frag[:i + 1]):
                    intra_molecular = True

            if intra_molecular:
                continue

            c6 = line.split()[2]

            # convert c6 constant from Haretrees/Bohr^6 to kcal/mol/A^6
            c6 *= constants.au_per_bohr6_to_kcal_per_ang6

            c6_constants.append(c6)

        # TODO: AVERAGE EQUIVELENT C6 CONSTANTS

        # read lines until we read the line before where charges are specified
        while True:
            try:
                qchem_out.readline().index("Ground-State ChElPG Net Atomic Charges")
                break
            except IndexError:
                qchem_out.readline() for i in range(3) # skip 3 more lines
                pass

        charges = []

        # read the charges for each atom
        for i in range(number_of_atoms):
            charge = float(qchem_out.readline().split()[2])
            charges.append(charge)

        return charges, polarizabilities, c6_constants
