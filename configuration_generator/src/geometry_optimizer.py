import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../qm_mb_energy_calculator/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import qcalc, output_writer
import settings_reader
import molecule_parser

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

    settings = settings_reader.SettingsReader(settings_path)

    # read the unoptimized geometry
    unopt_molecule = molecule_parser.xyz_to_molecules(unopt_path, settings)[0]

    # optimize the geometry
    opt_molecule, energy = qcalc.optimize(settings, unopt_molecule, settings.get("config_generator", "method"), settings.get("config_generator", "basis"))

    # write the optimized geometry to the output file
    with open(opt_path, "w") as opt_file:
        opt_file.write("{}\n".format(opt_molecule.get_num_atoms()))
        opt_file.write("{}\n".format(energy))
        opt_file.write("{}\n".format(opt_molecule.to_xyz()))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python geometry_optimizer.py <settings_path> <unopt_path> <opt_path>")
        exit(1)
    optimize_geometry(sys.argv[1], sys.argv[2], sys.argv[3])
