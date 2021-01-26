import itertools
from .psi4_calculator import Psi4Calculator
from .qchem_calculator import QchemCalculator
from mbfit.utils import SettingsReader, files, constants
from mbfit.exceptions import NoSuchLibraryError, LibraryCallError, PotentialFittingError
from mbfit.molecule import parse_training_set_file
from .model import Model

def get_calculator(settings_path, logging = True):
    """
    Gets a new Calculator object that can be used to perform QM calculations.

    Args:
        settings_path       - Local path to '.ini' settings file with all relevant settings.
        logging             - Not used at this time.

    Returns:
        A new Calculator object.
    """
    settings = SettingsReader(settings_path)
    if settings.get("energy_calculator", "code") == "psi4":
        return Psi4Calculator(settings_path, logging)
    elif settings.get("energy_calculator", "code") == "qchem":
        return QchemCalculator(settings_path, logging)
    else:
        raise NoSuchLibraryError(settings.get("energy_calculator", "code"))

def fill_energies(settings_path, input_configs_path, monomer_settings_paths, optimized_geometry_paths, output_configs_path, method, basis, cp):
    calculator = get_calculator(settings_path)

    files.init_file(output_configs_path)
    print("Calculating energies of optimized geometries...")
    counter = 0

    optimized_energies = []
    for monomer_settings_path, optimized_geometry_path in zip(monomer_settings_paths, optimized_geometry_paths):
        molecule = list(parse_training_set_file(optimized_geometry_path, SettingsReader(monomer_settings_path)))[0]
        try:
            model = Model(method, basis, False)

            # calculate the missing energy
            energy, log_path = calculator.calculate_energy(molecule, model, [0])
            optimized_energies.append(energy * constants.au_to_kcal)

        except LibraryCallError as e:
            raise PotentialFittingError("Energy calculation failed for optimized monomer {} with method {} and basis {}".format(molecule.get_name(), method, basis))
        counter += 1
        print("Completed optimized energy calculation for fragment number {}".format(counter))

    calc_successes, calc_failures = 0, 0
    total_successes, total_failures = 0, 0

    for molecule in parse_training_set_file(input_configs_path, SettingsReader(settings_path)):

        energies_dict = {}
        success = True

        for frag_indices, use_cp in get_energies_to_calculate(molecule.get_num_fragments(), cp):
            try:
                model = Model(method, basis, use_cp)

                # calculate the missing energy
                energy, log_path = calculator.calculate_energy(molecule, model, frag_indices)
                calc_successes += 1
                energies_dict[(frag_indices, use_cp)] = energy * constants.au_to_kcal

            except LibraryCallError as e:
                calc_failures += 1
                success = False

        if not success:
            total_failures += 1
            continue

        total_successes += 1

        deformation_energies = []

        for optimized_energy, deformed_energy in zip(optimized_energies, [item[1] for item in energies_dict.items() if len(item[0][0]) == 1 and item[0][1] is False]):
            deformation_energies.append(deformed_energy - optimized_energy)

        if molecule.get_num_fragments() == 1:
            interaction_energy = deformation_energies[0]
            binding_energy = deformation_energies[0]
        else:

            interaction_energy = energies_dict[(tuple(range(molecule.get_num_fragments())), False)]

            for m, energy in [(len(item[0][0]), item[1]) for item in energies_dict.items() if len(item[0][0]) < molecule.get_num_fragments() and item[0][1] is cp]:
                if (molecule.get_num_fragments() - m) % 2 == 1:
                    interaction_energy -= energy
                else:
                    interaction_energy += energy

            binding_energy = interaction_energy

            for deformed_energy in deformation_energies:
                binding_energy += deformed_energy

        with open(output_configs_path, "a") as output_configs_file:
            output_configs_file.write("{}\n".format(molecule.get_num_atoms()))
            output_configs_file.write("{} {}\n".format(binding_energy, interaction_energy))
            output_configs_file.write(molecule.to_xyz())
            output_configs_file.write("\n")

        if (total_successes + total_failures) % 10 == 0:
            print("{} Geometries complete!".format(total_successes + total_failures))

    print("Completed finding energies in training set. {} configurations included in training set, {} configurations not included".format(total_successes, total_failures)
          + " due to at least one failed calculation.")


def get_energies_to_calculate(num_bodies, cp):
    permutations = []
    for i in range(1, num_bodies + 1):
        perms = list(itertools.combinations(range(num_bodies), i))
        for p in perms:
            if cp and i < num_bodies:
                yield p, True
            yield p, False
