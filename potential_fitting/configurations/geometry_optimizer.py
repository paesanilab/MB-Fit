# absolute module imports
from potential_fitting import calculator
from potential_fitting.utils import SettingsReader
from potential_fitting.molecule import xyz_to_molecules

def optimize_geometry(settings_path, unopt_path, opt_path):
    """
    Optimizes a given geometry.

    Args:
        settings_path       - Local path to the ".ini" file with all relevent settings.
        unopt_path          - Local path to the ".xyz" file to read the unoptimized geometry from.
        opt_path            - Local path to the ".xyz" file to write the optimized geometry to.

    Returns:
        None.
    """

    settings = SettingsReader(settings_path)

    # read the unoptimized geometry
    unopt_molecule = xyz_to_molecules(unopt_path, settings)[0]

    # optimize the geometry
    opt_molecule, energy = calculator.optimize(settings, unopt_molecule, settings.get("config_generator", "method"), settings.get("config_generator", "basis"))

    # write the optimized geometry to the output file
    with open(opt_path, "w") as opt_file:
        opt_file.write("{}\n".format(opt_molecule.get_num_atoms()))
        opt_file.write("{}\n".format(energy))
        opt_file.write("{}\n".format(opt_molecule.to_xyz()))
