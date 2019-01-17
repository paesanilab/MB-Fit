# absolute package imports
import math, sqlite3

# absolute module imports
from potential_fitting.utils import constants
from potential_fitting.exceptions import NoEnergiesError, NoOptimizedEnergyError, MultipleOptimizedEnergiesError, NoEnergyInRangeError

# local module imports
from .database import Database

def generate_1b_training_set(settings_file, database_name, training_set_path, molecule_name, method, basis, cp, tag, e_min = 0, e_max = float('inf')):
    """
    Creates a 1b training set file from the calculated energies in a database.

    Args:
        settings_file       - Local path to the ".ini" file with all relevent settings information.
        database_name       - Local path to the database file. ".db" will be added if it does not already end in ".db".
        training_set_path   - Local path to file to write training set to.
        molecule_name       - The name of the molecule to generate a training set for.
        method              - Use energies calculated with this method. Use % for any method.
        basis               - Use energies calculated with this basis. Use % for any basis.
        cp                  - Use energies calculated with this cp. Use 0 for False, 1 for True, or % for any cp.
        tag                 - Use energies marked with this tag. Use % for any tag.
        e_min               - The minimum (inclusive) energy of any configuration to include in the training set.
        e_max               - The maximum (exclusive) energy of any configuration to include in the training set.

    Return:
        None.
    """
    
    # open the database
    with Database(database_name) as database:

        print("Creating a fitting input file from database {} into file {}".format(database_name, training_set_path))

        #intializing a counter
        count_configs = 0

        # get list of all [molecule, energies] pairs calculated in the database
        molecule_energy_pairs = list(database.get_energies(molecule_name, method, basis, cp, tag))

        # if there are no calculated energies, error and exit
        if len(molecule_energy_pairs) == 0:
            raise NoEnergiesError(database.file_name, molecule_name, method, basis, cp, tag)
        
        # find the optimized geometry energy from the database
        try:
            all_opt_energies = list(database.get_energies(molecule_name, method, basis, cp, tag, True))
            if len(all_opt_energies) > 1:
                raise MultipleOptimizedEnergiesError(database.file_name, molecule_name, method, basis, cp, tag, 
                        len(all_opt_energies))
            opt_energies = all_opt_energies[0][1]
        except IndexError:
            raise NoOptimizedEnergyError(database.file_name, molecule_name, method, basis, cp, tag) from None

        # open output file for writing
        with open(training_set_path, "w") as output:

            # loop thru each molecule, energy pair
            for molecule_energy_pair in molecule_energy_pairs:
                molecule = molecule_energy_pair[0]
                energies = molecule_energy_pair[1]

                # monomer deformation energy
                deformation_energy_hartrees = (float(energies[0]) - float(opt_energies[0]))
                deformation_energy_kcalmol = deformation_energy_hartrees * constants.au_to_kcal # Converts Hartree to kcal/mol
                if deformation_energy_kcalmol < e_max and deformation_energy_kcalmol - e_min > -0.000000000001:

                    # write the number of atoms to the output file
                    output.write(str(molecule.get_num_atoms()) + "\n")
                    output.write(str(deformation_energy_kcalmol) + " ") 
                    output.write("\n")

                    # write the molecule's atoms' coordinates to the xyz file
                    output.write(molecule.to_xyz() + "\n")

                    #increment the counter
                    count_configs += 1

            if count_configs == 0:
                raise NoEnergyInRangeError(database.file_name, molecule_name, method, basis, cp, tag, e_min, e_max)
            print("Generated training set with " + str(count_configs) + " Configurations.")



                    

def generate_2b_training_set(settings, database_name, training_set_path, monomer_1_name, monomer_2_name, method, basis,
        cp, tag, e_bind_max = float('inf'), e_mon_max = float('inf')):
    """"
    Creates a 2b training set file from the calculated energies in a database.

    Args:
        settings_file       - Local path to the ".ini" file with all relevent settings information.
        database_name       - Local path to the database file. ".db" will be added if it does not already end in ".db".
        training_set_path   - Local path to file to write training set to.
        monomer_1_name      - The name of the first monomer in the dimer.
        monomer_2_name      - The name of the second monomer in the dimer.
        method              - Use energies calculated with this method. Use % for any method.
        basis               - Use energies calculated with this basis. Use % for any basis.
        cp                  - Use energies calculated with this cp. Use 0 for False, 1 for True, or % for any cp.
        tag                 - Use energies marked with this tag. Use % for any tag.
        e_bind_max          - Maximum binding energy allowed
        e_mon_max           - Maximum monomer deformation energy allowed

    Return:
        None.
    """
    
    # open the database
    with Database(database_name) as database:

        print("Creating a fitting input file from database {} into file {}".format(database_name, training_set_path))

        # construct name of molecule from name of monomers
        molecule_name = "-".join([monomer_1_name, monomer_2_name])
        # get list of all [molecule, energies] pairs calculated in the database
        molecule_energy_pairs = list(database.get_energies(molecule_name, method, basis, cp, tag))

        # if there are no calculated energies, error and exit
        if len(molecule_energy_pairs) == 0:
            raise NoEnergiesError(database.file_name, molecule_name, method, basis, cp, tag)
        
        # find the optimized geometry energy of the two monomers from the database
        try:
            all_opt_energies = list(database.get_energies(monomer_1_name, method, basis, cp, tag, True))
            if len(all_opt_energies) > 1:
                raise MultipleOptimizedEnergiesError(database.file_name, monomer_1_name, method, basis, cp, tag, 
                        len(all_opt_energies))
            monomer_1_opt_energy = all_opt_energies[0][1][0]
        except IndexError:
            raise NoOptimizedEnergyError(database.file_name, monomer_1_name, method, basis, cp, tag)

        try:
            all_opt_energies = list(database.get_energies(monomer_2_name, method, basis, cp, tag, True))
            if len(all_opt_energies) > 1:
                raise MultipleOptimizedEnergiesError(database.file_name, monomer_2_name, method, basis, cp, tag,
                        len(all_opt_energies))
            monomer_2_opt_energy = all_opt_energies[0][1][0]
        except IndexError:
            raise NoOptimizedEnergyError(database.file_name, monomer_2_name, method, basis, cp, tag)

        # open output file for writing
        with open(training_set_path, "w") as output:

            for molecule_energy_pair in molecule_energy_pairs:
                molecule = molecule_energy_pair[0]
                energies = molecule_energy_pair[1]

                
                # calculate the interaction energy of the dimer as E01 - E0 - E1 all computed in the dimer basis set
                # if cp is set to true (otherwise in their own basis set)
                interaction_energy = (energies[2] - energies[1] - energies[0]) * constants.au_to_kcal

                # if there are more energies than for a non-cp enabled calculation, we know that this is a cp enabled 
                # alculation
                if len(energies) > 3:

                    # if cp is enabled for this calculation, then computed the monomer deformation energies as
                    # deformed energy - optimized energy (in the monomer basis set)
                    monomer1_energy_deformation = (energies[3] - monomer_1_opt_energy) * constants.au_to_kcal
                    monomer2_energy_deformation = (energies[4] - monomer_2_opt_energy) * constants.au_to_kcal

                # otherwise it is a non-cp enabled calculation
                else:

                    # if cp is not enabled for this calculation, then computed the monomer deformation energies as
                    # deformed energy - optimized energy (in the monomer basis set)
                    monomer1_energy_deformation = (energies[0] - monomer_1_opt_energy) * constants.au_to_kcal
                    monomer2_energy_deformation = (energies[1] - monomer_2_opt_energy) * constants.au_to_kcal

                # calculate binding energy
                binding_energy = interaction_energy - monomer1_energy_deformation - monomer2_energy_deformation

                # write the configuration only if is below the thresholds
                if (binding_energy < e_bind_max) and (monomer1_energy_deformation < e_mon_max) and (monomer2_energy_deformation < e_mon_max):
                    # write the number of atoms to the output file 
                    output.write(str(molecule.get_num_atoms()) + "\n")
    
                    output.write("{} {} {} {}".format(binding_energy, interaction_energy, monomer1_energy_deformation,
                            monomer2_energy_deformation))
                
                    output.write("\n")
    
                    # write the molecule's atoms' coordinates to the xyz file
                    output.write(molecule.to_xyz() + "\n")
