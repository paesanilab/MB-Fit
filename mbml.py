import sys, os
import potential_fitting
from potential_fitting import calculator
from potential_fitting.utils import SettingsReader

# check arguments!

if len(sys.argv) != 3:
    print("Incorrect arguments!")
    print("Usage: python3 mbml <training_set_path> <settings_path>")
    sys.exit(1)

training_set_path = sys.argv[1]
settings_path = sys.argv[2]



settings = SettingsReader(settings_file)

# STEP 1: calculate energies in the training set.

calculator.fill_training_set(settings_path, training_set_path, settings.get("files", "filled_energies_path")
                                                             , settings.get("model", "method")
                                                             , settings.get("model", "basis")
                                                             , settings.get("model", "cp"))