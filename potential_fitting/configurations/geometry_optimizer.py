# absolute module imports
from potential_fitting import calculator
from potential_fitting.calculator import Model
from potential_fitting.utils import SettingsReader, files, system
from potential_fitting.molecule import xyz_to_molecules

def optimize_geometry(settings_path, unopt_path, opt_path, method, basis, arguments={}):
    """
    Optimizes a given geometry.

    Args:
        settings_path       - Local path to the ".ini" file with all relevent settings.
        unopt_path          - Local path to the ".xyz" file to read the unoptimized geometry from.
        opt_path            - Local path to the ".xyz" file to write the optimized geometry to.
        arguments           - Dictionary of extra arguments to be passed to the QM code doing the calculation.

    Returns:
        None.
    """
    system.format_print("Beginning geometry optimization of {}...".format(unopt_path),
            bold=True, color=system.Color.YELLOW)

    settings = SettingsReader(settings_path)

    calc = calculator.get_calculator(settings_path)
    model = Model(method, basis)

    # read the unoptimized geometry
    unopt_molecule = xyz_to_molecules(unopt_path, settings)[0]

    # optimize the geometry
    opt_molecule, energy, log_file = calc.optimize_geometry(unopt_molecule, model, arguments=arguments)

    opt_path = files.init_file(opt_path, files.OverwriteMethod.get_from_settings(settings))

    # write the optimized geometry to the output file
    with open(opt_path, "w") as opt_file:
        opt_file.write("{}\n".format(opt_molecule.get_num_atoms()))
        opt_file.write("{}\n".format(energy))
        opt_file.write("{}\n".format(opt_molecule.to_xyz()))

    system.format_print("Geometry optimization complete! Optimized geometry in {}.".format(opt_path),
            bold=True, color=system.Color.GREEN)
