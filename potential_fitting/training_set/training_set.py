from potential_fitting.exceptions import XYZFormatError
from potential_fitting.molecule import parse_training_set_file

from . import TrainingSetElement

class TrainingSet:

    @staticmethod
    def get_training_set_from_data(molecules, **energies_dict):
        elements = []

        for index, molecule in enumerate(molecules):

            molecule_energies = {}

            for key, energies in energies_dict.items():
                molecule_energies[key] = energies[index]

            element = TrainingSetElement(molecule, **molecule_energies)
            elements.append(element)

        return TrainingSet(elements)

    @staticmethod
    def get_training_set_from_xyz_file(path_to_xyz_file, settings, energy_names, is_training_format = True):

        molecules = list(parse_training_set_file(path_to_xyz_file, settings=settings))

        num_atoms = sum([int(atom_count) for atom_count in settings.get("molecule", "fragments").split(",")])

        energies_dict = {}

        for energy_name in energy_names:
            energies_dict[energy_name] = []

        with open(path_to_xyz_file, "r") as xyz_file:
            while True:
                atom_count_line = xyz_file.readline()
                while len(atom_count_line.split()) == 0:
                    if atom_count_line == "":
                        return TrainingSet.get_training_set_from_data(molecules, **energies_dict)

                    atom_count_line = xyz_file.readline()

                energies_line = xyz_file.readline()
                if energies_line == "":
                    raise XYZFormatError("ran out of lines to read from xyz file {} in the middle of a molecule".format(xyz_file.name), "make sure atoms_per_fragment, the atom count line in your xyz file, and the number of atom lines in your xyz file all agree.")

                energies = [float(e) for e in energies_line.split()]

                if len(energies) != len(energy_names):
                    if is_training_format:
                        raise XYZFormatError("Expected {} enegies in line '{}' but found only {}.".format(len(energy_names), energies_line, len(energies)))
                    else:
                        energies = [0.0]*len(energy_names)

                for index, energy_name in enumerate(energy_names):
                    energies_dict[energy_name].append(energies[index])

                for index in range(num_atoms):
                    xyz_file.readline()

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

    def add_energies(self, energy_name, energies):

        for element, energy in zip(self.get_elements(), energies):
            element.add_energy(energy_name, energy)

    def __str__(self):
        return "\n".join(str(element) for element in self.get_elements())
