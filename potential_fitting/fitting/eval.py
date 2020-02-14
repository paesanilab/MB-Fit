import os

from potential_fitting.utils import system, files
from potential_fitting.training_set import TrainingSet

class Evaluator:

    def __init__(self, settings, path_to_eval_file):
        self.settings = settings
        self.path_to_eval_file = path_to_eval_file

    def eval(self, path_to_training_set_file, settings):

        correlation_file_path = files.init_file(self.settings.get("files", "logs"), os.path.join("eval_correlation.dat"))

        with open(correlation_file_path, "w") as correlation_file:
            system.call(self.path_to_eval_file, path_to_training_set_file, out_file=correlation_file)

        with open(correlation_file_path, "r") as correlation_file:

            fit_energies = []
            for line in correlation_file.readlines():
                if not line.startswith("#"):
                    fit_energies.append(float(line.split()[1]))

        training_set = TrainingSet.get_training_set_from_xyz_file(path_to_training_set_file, settings, energy_names=["weight_energy", "ref_energy"])

        training_set.add_energies("fit_energies", fit_energies)

        return training_set