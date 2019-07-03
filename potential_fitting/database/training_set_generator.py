# absolute module imports
from potential_fitting.utils import constants, SettingsReader, files
from potential_fitting.exceptions import NoEnergiesError, NoOptimizedEnergyError, MultipleOptimizedEnergiesError, NoEnergyInRangeError

# local module imports
from .database import Database


def generate_1b_training_set(settings_path, database_config_path, training_set_path, molecule_name, method, basis, cp, *tags, e_min=0, e_max=float('inf')):
    """
    Writes a 1b training set to the given file from the calculated energies in a database.

    Args:
        settings_path       - Local path to the ".ini" file with all relevent settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        training_set_path   - Local path to file to write training set to.
        molecule_name       - The name of the molecule to generate a training set for.
        method              - Use energies calculated with this method. Use % for any method.
        basis               - Use energies calculated with this basis. Use % for any basis.
        cp                  - Use energies calculated with this cp. Use 0 for False, 1 for True, or % for any cp.
        tags                - Use energies marked with one or more of these tags. Use % for any tag.
        e_min               - The minimum (inclusive) energy of any configuration to include in the training set.
        e_max               - The maximum (exclusive) energy of any configuration to include in the training set.

    Return:
        None.
    """

    settings = SettingsReader(settings_path)

    names = settings.get("molecule", "names").split(",")
    SMILES = settings.get("molecule", "SMILES").split(",")
    
    # open the database
    with Database(database_config_path) as database:

        print("Creating a fitting input file from database into file {}".format(training_set_path))

        #intializing a counter
        count_configs = 0

        with open(files.init_file(training_set_path, files.OverwriteMethod.get_from_settings(settings)), "w") as output:

            for molecule, deformation_energy_hartrees in database.get_1B_training_set(molecule_name, names, SMILES, method, basis, cp, *tags):

                deformation_energy_kcalmol = deformation_energy_hartrees * constants.au_to_kcal # Converts Hartree to kcal/mol
                if deformation_energy_kcalmol < e_max and deformation_energy_kcalmol - e_min > -0.000000000001:
                    # write the number of atoms to the output file
                    output.write(str(molecule.get_num_atoms()) + "\n")
                    output.write(str(deformation_energy_kcalmol) + " ")
                    output.write("\n")

                    # write the molecule's atoms' coordinates to the xyz file
                    output.write(molecule.to_xyz() + "\n")

                    # increment the counter
                    count_configs += 1

            if count_configs == 0:
                raise Exception

            print("Generated training set with " + str(count_configs) + " Configurations.")


def generate_2b_training_set(settings_path, database_config_path, training_set_path, molecule_name, method, basis,
        cp, *tags, e_bind_max=float('inf'), e_mon_max=float('inf')):
    """"
    Creates a 2b training set file from the calculated energies in a database.

    Args:
        settings_path       - Local path to the ".ini" file with all relevent settings information.
        database_config_path - .ini file containing host, port, database, username, and password.
                    Make sure only you have access to this file or your password will be compromised!
        training_set_path   - Local path to file to write training set to.
        molecule_name       - The name of this dimer.
        monomer1_name       - The name of one monomer in the dimer.
        monomer2_name       - The name of the other monomer in the dimer.
        method              - Use energies calculated with this method. Use % for any method.
        basis               - Use energies calculated with this basis. Use % for any basis.
        cp                  - Use energies calculated with this cp. Use 0 for False, 1 for True, or % for any cp.
        tags                - Use energies marked with at least one of these tags. Use % for any tag.
        e_bind_max          - Maximum binding energy allowed
        e_mon_max           - Maximum monomer deformation energy allowed

    Return:
        None.
    """

    settings = SettingsReader(settings_path)

    names = settings.get("molecule", "names").split(",")
    SMILES = settings.get("molecule", "SMILES").split(",")
    
    # open the database
    with Database(database_config_path) as database:

        print("Creating a fitting input file from database into file {}".format(training_set_path))

        # initializing a counter
        count_configs = 0

        with open(files.init_file(training_set_path, files.OverwriteMethod.get_from_settings(settings)), "w") as output:
            for molecule, binding_energy, interaction_energy, monomer1_energy_deformation, monomer2_energy_deformation in database.get_2B_training_set(molecule_name, names, SMILES, method, basis, cp, *tags):

                binding_energy *= constants.au_to_kcal
                interaction_energy *= constants.au_to_kcal
                monomer1_energy_deformation *= constants.au_to_kcal
                monomer2_energy_deformation *= constants.au_to_kcal

                # write the configuration only if is below the thresholds
                if (binding_energy < e_bind_max) and (monomer1_energy_deformation < e_mon_max) and (
                        monomer2_energy_deformation < e_mon_max):

                    # write the number of atoms to the output file
                    output.write(str(molecule.get_num_atoms()) + "\n")

                    output.write("{} {} {} {}".format(binding_energy, interaction_energy, monomer1_energy_deformation,
                                                      monomer2_energy_deformation))

                    output.write("\n")

                    # write the molecule's atoms' coordinates to the xyz file
                    output.write(molecule.to_xyz() + "\n")

                    # increment the counter
                    count_configs += 1

        if count_configs == 0:
            raise Exception

        print("Generated training set with " + str(count_configs) + " Configurations.")
