from potential_fitting import calculator
from potential_fitting.utils import SettingsReader
from potential_fitting.molecule import xyz_to_molecules

def optimize_geometry(settings_path, unopt_path, opt_path):
    """
    Optimizes a given geometry

    Args:
        settings_path - path to the .ini file with all relevent settings
        unopt_path  - path to the geometry .xy file to optimize
        opt_path    - path to .xyz file to write the optimized geometry

    Returns:
        None
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
