import sys, os
import potential_fitting
from potential_fitting import calculator
from potential_fitting.utils import SettingsReader

# check arguments!

if len(sys.argv) != 3:
    print("Incorrect arguments!")
    print("Usage: python3 mbml <training_configs_path> <settings_path>")
    sys.exit(1)

input_set_path = sys.argv[1]
settings_path = sys.argv[2]

settings = SettingsReader(settings_path)

monomer_settings_paths = settings.getlist("files", "monomer_settings_paths")


names = settings.get("molecule", "names").split(",")
fragments = settings.get("molecule", "fragments").split(",")
charges = settings.get("molecule", "charges").split(",")
spins = settings.get("molecule", "spins").split(",")
symmetries = settings.get("molecule", "symmetry").split(",")
SMILES = settings.get("molecule", "SMILES").split(",")

for monomer_settings_path, name, fragment, charge, spin, symmetry, SMILE in zip(monomer_settings_paths, names, fragments, charges, spins, symmetries, SMILES):
    monomer_setting = SettingsReader(settings_path)
    monomer_setting.set("molecule", "names", name)
    monomer_setting.set("molecule", "fragments", fragment)
    monomer_setting.set("molecule", "charges", charge)
    monomer_setting.set("molecule", "spins", spin)
    monomer_setting.set("molecule", "symmetry", symmetry)
    monomer_setting.set("molecule", "SMILES", SMILE)
    monomer_setting.write(monomer_settings_path)

optimized_geometry_paths = settings.getlist("files", "optimized_geometry_paths")
training_set_path = settings.get("files", "training_set_path")

# STEP 1: calculate energies in the training set.
calculator.fill_energies(settings_path, input_set_path, monomer_settings_paths, optimized_geometry_paths, training_set_path,
                         settings.get("model", "method"),
                         settings.get("model", "basis"),
                         settings.get("model", "cp"))