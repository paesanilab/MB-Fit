import itertools
from .psi4_calculator import Psi4Calculator
from .qchem_calculator import QchemCalculator
from potential_fitting.utils import SettingsReader, files
from potential_fitting.exceptions import NoSuchLibraryError, LibraryCallError
from potential_fitting.molecule import parse_training_set_file
from .model import Model
def get_calculator(settings_path, logging = True):
    settings = SettingsReader(settings_path)
    if settings.get("energy_calculator", "code") == "psi4":
        return Psi4Calculator(settings_path, logging)
    elif settings.get("energy_calculator", "code") == "qchem":
        return QchemCalculator(settings_path, logging)
    else:
        raise NoSuchLibraryError(settings.get("energy_calculator", "code"))

def fill_energies(settings_path, input_configs_path, output_configs_path, method, basis, cp):
    calculator = get_calculator(settings_path)

    files.init_file(output_configs_path)

    for molecule in parse_training_set_file(input_configs_path, settings_path):

        energies_dict = {}

        for frag_indices, use_cp in get_energies_to_calculate(molecule.get_fragments(), cp):
            try:
                model = Model(method, basis, use_cp)

                # calculate the missing energy
                energy, log_path = calculator.calculate_energy(molecule, model, frag_indices)
                successes += 1

            except LibraryCallError as e:
                failures += 1

def get_energies_to_calculate(num_bodies, cp):
    permutations = []
    for i in range(1, num_bodies + 1):
        perms = list(itertools.combinations(range(num_bodies), i))
        for p in perms:
            if cp and i < num_bodies:
                yield p, True
            yield p, False



