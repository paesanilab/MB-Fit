from . import TrainingSetElement

class TrainingSet:

    @staticmethod
    def get_training_set_from_data(molecules, **energies):
        elements = []

        for index, molecule in enumerate(molecules):

            molecule_energies = {}

            for key, energies in energies:
                molecule_energies[key] = energies[index]

            element = TrainingSetElement(molecule, **molecule_energies)
            elements.append(element)

        return TrainingSet(elements)

    def __init__(self, elements):
        self.elements = elements

    def get_elements(self):
        return self.elements

    def get_molecules(self):
        mols = []
        for element in self.get_elements():
            mols.append(element.get_molecule())

        return mols

    def has_energies(self, energy_name):
        for element in self.get_elements():
            if not element.has_energy(energy_name):
                return False

        return True

    def get_energies(self, energy_name):
        energies = []
        for element in self.get_elements():
            energies.append(element.get_energy(energy_name))

        return energies

    def split_at_threshold(self, energy_name, threshold):

        low_elements = []
        high_elements = []

        for element in self.get_elements():
            if element.get_energy(energy_name) < threshold:
                low_elements.append(element)
            else:
                high_elements.append(element)

        return TrainingSet(low_elements), TrainingSet(high_elements)
