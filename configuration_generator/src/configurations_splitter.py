import sys, os
import random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../qm_mb_energy_calculator/src")
import molecule_parser
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import settings_reader

def splitConfigurations(settings_file, configurations_path, training_set_path, test_set_path, training_set_size, molecular_descriptor = None):
    """
    Splits a set of configurations into a training set and a test set using furthest
    point sampling by some measure defined by the given molecular_descriptor

    Args:
        settings_file - path to .ini file containing relevent settings information
        configurations_path - file path to the xyz file with configurations to be divided into a training and test set
        training_set_path - file path to the xyz file to write the training set to.
        test_set_path - file path to the xyz file to write the test set to.
        training_set_size - the desired size of the training set, all other molecules will be put into the test set
        molecular_descriptor - the MolecularDescriptor used to measure the distance between two molecules.
    """
    
    if molecular_descriptor is None:
        molecular_descriptor = RMSDDescriptor()
    molecules = molecule_parser.xyz_to_molecules(configurations_path, settings_reader.SettingsReader(settings_file))

    training_set = []

    if training_set_size > 0:
        first_train_set_molecule = random.choice(molecules)
        molecules.remove(first_train_set_molecule)
        training_set.append(first_train_set_molecule)

    for i in range(training_set_size - 1):
        minimum_difference_to_train_set = []

        for molecule in molecules:
            minimum_difference_to_train_set.append(molecular_descriptor.difference(molecule, training_set[0]))

        for molecule_index, molecule in enumerate(molecules):

            for train_set_molecule in training_set[1:]:

                difference = molecular_descriptor.difference(molecule, train_set_molecule)

                if difference < minimum_difference_to_train_set[molecule_index]:
                    minimum_difference_to_train_set[molecule_index] = difference

        largest_min_difference_molecule_index = 0

        for molecule_index in range(len(molecules)):
            if minimum_difference_to_train_set[molecule_index] > minimum_difference_to_train_set[largest_min_difference_molecule_index]:
                largest_min_difference_molecule_index = molecule_index

        train_set_molecule = molecules[largest_min_difference_molecule_index]
        training_set.append(train_set_molecule)
        molecules.remove(train_set_molecule)

    with open(training_set_path, "w") as training_set_file:

        for molecule in training_set:

            # write number of atoms to training set file
            training_set_file.write("{}\n".format(molecule.get_num_atoms()))

            # write comment line to training set file
            training_set_file.write("\n")

            # write geometry to training set file
            training_set_file.write("{}\n".format(molecule.to_xyz()))
        
    with open(test_set_path, "w") as test_set_file:

        for molecule in molecules:

            # write number of atoms to test set file
            test_set_file.write("{}\n".format(molecule.get_num_atoms()))

            # write comment line to test set file
            test_set_file.write("\n")

            # write geometry to test set file
            test_set_file.write("{}\n".format(molecule.to_xyz()))

class MolecularDescriptor():
    def difference(molecule1, molecule2):
        raise NotImplementedError

class RMSDDescriptor(MolecularDescriptor):
    def difference(self, molecule1, molecule2):
        return molecule1.rmsd(molecule2)
