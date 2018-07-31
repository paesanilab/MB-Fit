import molecule
import configparser

def optimize_geometry(settings, unopt_path, opt_path):

    config = configparser.ConfigParser(allow_no_value=False)
    config.read(settings)
    
    with open(unopt_path, "r") as unopt_file:
        unopt_molecule = molecule.Molecule(unopt_file.read())

    opt_molecule, energy = qcalc.optimize(unopt_molecule, config)

    output_writer.write_optimized_geo(opt_molecule, energy, opt_path)

