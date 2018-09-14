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

    All the molecules will be moved to their center of mass and rotated on their principal axes in the process

    Args:
        settings_file - path to .ini file containing relevent settings information
        configurations_path - file path to the xyz file with configurations to be divided into a training and test set
        training_set_path - file path to the xyz file to write the training set to.
        test_set_path - file path to the xyz file to write the test set to.
        training_set_size - the desired size of the training set, all other molecules will be put into the test set
        molecular_descriptor - the MolecularDescriptor used to measure the difference between two molecules.
    """

    # if the user did not specify a molecular descriptor, use the RMSD descriptor    
    if molecular_descriptor is None:
        molecular_descriptor = RMSDDescriptor()

    # parse the configurations into a list of [id, molecule]
    molecules = list(enumerate(molecule_parser.xyz_to_molecules(configurations_path, settings_reader.SettingsReader(settings_file))))

    # move all molecules to their center of mass and rotate them on their principal axes
    for index, molecule in molecules:
        molecule.move_to_center_of_mass()
        molecule.rotate_on_principal_axes()

    # construct a matrix where the item at index [i, k] is the difference between the ith molecule and the kth molecule
    # it starts empty, and will be filled in
    difference_matrix = [[None for k in range(i + 1)] for i in range(len(molecules) - 1)]

    # this array will hold the [id, molecule] pairs of the training set, those for the test set will be left in molecules
    training_set = []

    # if training_set_size is at least 1, then move a random molecule from molecules to training_set
    if training_set_size > 0:
        first_train_set_molecule = random.choice(molecules)
        molecules.remove(first_train_set_molecule)
        training_set.append(first_train_set_molecule)

    # loop once for every molecule to be added to the training_set
    for i in range(training_set_size - 1):

        # this array holds a list of the minimum differences to any element in the training set for each element in molecules
        minimum_difference_to_train_set = []

        # loop over each [id, molecule] pair still in molecules
        for index, molecule in molecules:

            # if the difference between this molecule and the first molecule in the training_set has not been calculated, calculate it.
            if difference_matrix[index - 1][0] is None:
                difference_matrix[index - 1][0] = molecular_descriptor.difference(molecule, training_set[0][1])

            # initialize this molecule's entry in minimum_difference_to_train_set to the difference between this molecule and the first molecule in the training_set
            minimum_difference_to_train_set.append(difference_matrix[index - 1][0])

        # loop over each [id, molecule] pair still in molecules
        for molecule_index, index_molecule_pair in enumerate(molecules):

            # loop over each [id, molecule] pair in the training_set
            for train_index_molecule_pair in training_set[1:]:

                # sort the two ids of the molecules being compared so the smaller one goes into index2 and the larger into index1
                index2, index1 = sorted([index_molecule_pair[0], train_index_molecule_pair[0]])
                index1 -= 1

                # if the difference for this pair has not been computed, compute it
                if difference_matrix[index1][index2] is None:
                    difference_matrix[index1][index2] = molecular_descriptor.difference(index_molecule_pair[1], train_index_molecule_pair[1])

                # loop up the difference between the two molecules in the difference_matrix
                difference = difference_matrix[index1][index2]

                # if this difference is less than the current minimum difference for this molecule, update the minimum difference
                if difference < minimum_difference_to_train_set[molecule_index]:
                    minimum_difference_to_train_set[molecule_index] = difference

        # this holds the index of the molecule with the largest minimum difference to the training_set
        largest_min_difference_molecule_index = 0

        # loop over each molecule not yet in the training_set
        for molecule_index in range(len(molecules)):

            # if its minimum difference to the training_set is greater than the current furthest molecule's minimum difference to the training set, then this molecule becoms the new furthest molecule
            if minimum_difference_to_train_set[molecule_index] > minimum_difference_to_train_set[largest_min_difference_molecule_index]:
                largest_min_difference_molecule_index = molecule_index

        # move the molecule furthest from the training set to the training set
        train_set_molecule = molecules[largest_min_difference_molecule_index]
        training_set.append(train_set_molecule)
        molecules.remove(train_set_molecule)

    with open(training_set_path, "w") as training_set_file:

        for index, molecule in training_set:

            # write number of atoms to training set file
            training_set_file.write("{}\n".format(molecule.get_num_atoms()))

            # write comment line to training set file
            training_set_file.write("\n")

            # write geometry to training set file
            training_set_file.write("{}\n".format(molecule.to_xyz()))
        
    with open(test_set_path, "w") as test_set_file:

        for index, molecule in molecules:

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

class RMSDDistanceDescriptor(MolecularDescriptor):
    def difference(self, molecule1, molecule2):
        return molecule1.distancermsd(molecule2)

class RandomDescriptor(MolecularDescriptor):
    def difference(self, molecule1, molecule2):
        return random.random()
