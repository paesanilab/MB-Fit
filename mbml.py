import sys, os, argparse
import potential_fitting
from potential_fitting import calculator
from potential_fitting.utils import SettingsReader, files
from potential_fitting.molecule import Molecule

# check arguments!

parser = argparse.ArgumentParser(description="Some Description of the Code",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# file path arguments.

parser.add_argument('--opt_geo', '-g', dest='optimized_geometry_path', type=str, required=True,
                    help='File path to the ".xyz" file containing the optimized geometry for a monomer or the '
                         'optimized geometries for a 2b+ molecule.')
parser.add_argument('--training_set_input', '-ti', dest='training_set_input_path', type=str, required=True,
                    help='File path to the ".xyz" file containing the geometries to include in the training set.')
skip_training_set_calculations_group = parser.add_mutually_exclusive_group(required=True)
skip_training_set_calculations_group.add_argument('--training_set_output', '-to', dest='training_set_output_path',
                                                  type=str, required=False,
                                                  help='File path to where the ".xyz" file containing the geometries'
                                                       'and energies of the training set should be written.')
parser.add_argument('--log_files', '-l', dest='log_path', type=str, required=True,
                    help='File path to directory in which to store log files for all calculations.')
parser.add_argument('--properties_path', '-p', dest='properties_path',
                                                  type=str, required=True,
                                                  help='File path to where the ".ini" file containing the properties'
                                                       'of the optimized geometry is to be written.')
parser.add_argument('--1b_config_paths', '-c1', nargs='+', dest='config_1b_paths', type=str, required=False,
                    help='List of file paths to ".ini" files containing 1b properties for each monomer. Only required '
                         'if your molecule has 2 or more fragments')
parser.add_argument('--2b_config_paths', '-c2', nargs='+', dest='config_2b_paths', type=str, required=False,
                    help='List of file paths to ".ini" files containing 2b properties for each dimer. Only required '
                         'if your molecule has 3 or more fragments')

# Properties that define the molecule

parser.add_argument('--num_fragments', '-nf', dest='number_of_fragments', type=int, required=True,
                    help='The number of fragments in this molecule.')
parser.add_argument('--fragment_names', '-fn', dest='fragment_names', type=str, required=True,
                    help='Comma delimited list of the names of each fragment.')
parser.add_argument('--fragment_atoms', '-fa', dest='fragment_atoms', type=str, required=True,
                    help='Comma delimited list of the number of atoms in each fragment.')
parser.add_argument('--fragment_charges', '-fc', dest='fragment_charges', type=str, required=True,
                    help='Comma delimited list of the net charge of each fragment.')
parser.add_argument('--fragment_spin_multiplicities', '-fs', dest='fragment_spin_multiplicities', type=str, required=True,
                    help='Comma delimited list of the spin multiplicity of each fragment.')
parser.add_argument('--fragment_symmetries', '-fm', dest='fragment_symmetries', type=str, required=True,
                    help='Comma delimited list of the symmetry string of each fragment. For example "A1B2X2" or "A3B2".')
parser.add_argument('--fragment_smiles', '-fl', dest='fragment_smiles', type=str, required=True,
                    help='Comma delimited list of the SMILE string of each fragment. Must explicitly include each H.')

# The QM model

parser.add_argument('--method', '-mm', dest='model_method', type=str, required=True,
                    help='QM method to use for training set calculations.')
parser.add_argument('--basis', '-mb', dest='model_basis', type=str, required=True,
                    help='QM basis set to use for training set calculations.')
parser.add_argument('--counterpoise', '-mc', dest='model_counterpoise_correction', type=bool, required=True,
                    help='If True, the counterpoise correction will be used for training set calculations.')

# Library options

parser.add_argument('--code', '-c', dest='code', type=str, required=True,
                    help='What code to use for calculations, options are "psi4" or "qchem".')
parser.add_argument('--num_threads', '-nt', dest='num_threads', type=str, required=False, default=1,
                    help='Number of threads to use for operations that support multithreading.')

# Optional arguments to skip various parts of the procedure.

skip_training_set_calculations_group.add_argument('--skip_training_set_calculations', '-sc', dest='calculate_training_set', required=False, action='store_false', default=True,
                    help='If included, then the input training set will be assumed to include calculated energies already, so no'
                         'calculations will be performed to put energies into the training set.')
parser.add_argument('--skip_properties_calculations', '-sp', dest='calculate_properties', required=False, action='store_false', default=True,
                    help='If included, then no calculation will be run to find the optimized geometry\'s properties. Instead, the properties path'
                         'specified with "--properties_path" will be assumed to already contain the properties.')

args = parser.parse_args()

if args.config_1b_paths is None:
    args.config_1b_paths = []
if args.config_2b_paths is None:
    args.config_2b_paths = []

# construct the settings reader to write the settings.ini

settings = SettingsReader()

settings.set('files', 'log_path', args.log_path)

settings.set('molecule', 'names', args.fragment_names)
settings.set('molecule', 'fragments', args.fragment_atoms)
settings.set('molecule', 'charges', args.fragment_charges)
settings.set('molecule', 'spins', args.fragment_spin_multiplicities)
settings.set('molecule', 'symmetry', args.fragment_symmetries)
settings.set('molecule', 'SMILES', args.fragment_smiles)

settings.set('model', 'method', args.model_method)
settings.set('model', 'basis', args.model_basis)
settings.set('model', 'cp', args.model_counterpoise_correction)

settings.set('energy_calculator', 'code', args.code)
settings.set('qchem', 'num_threads', args.num_threads)
settings.set('psi4', 'num_threads', args.num_threads)

temp_file_path = os.path.join(settings.get('files', 'log_path'), "temp_files")
settings_file_path = os.path.join(temp_file_path, "settings.ini")

settings.write(settings_file_path)

monomer_settings_paths = [os.path.join(temp_file_path, "monomer{}.ini".format(i)) for i in range(args.number_of_fragments)]


names = settings.get("molecule", "names").split(",")
fragments = settings.get("molecule", "fragments").split(",")
charges = settings.get("molecule", "charges").split(",")
spins = settings.get("molecule", "spins").split(",")
symmetries = settings.get("molecule", "symmetry").split(",")
SMILES = settings.get("molecule", "SMILES").split(",")

for monomer_settings_path, name, fragment, charge, spin, symmetry, SMILE in zip(monomer_settings_paths, names, fragments, charges, spins, symmetries, SMILES):
    monomer_setting = SettingsReader(settings_file_path)
    monomer_setting.set("molecule", "names", name)
    monomer_setting.set("molecule", "fragments", fragment)
    monomer_setting.set("molecule", "charges", charge)
    monomer_setting.set("molecule", "spins", spin)
    monomer_setting.set("molecule", "symmetry", symmetry)
    monomer_setting.set("molecule", "SMILES", SMILE)
    monomer_setting.write(monomer_settings_path)

optimized_geometry_paths = [os.path.join(temp_file_path, "opt_geo{}.ini".format(i)) for i in range(args.number_of_fragments)]

opt_molecule = Molecule.read_xyz_path_direct(args.optimized_geometry_path, settings)

for i, optimized_geometry_path in enumerate(optimized_geometry_paths):
    with open(optimized_geometry_path, "w") as optimized_geometry_file:
        optimized_geometry_file.write("{}\n".format(opt_molecule.get_fragments()[0].get_num_atoms()))
        optimized_geometry_file.write("optimized geometry for fragment {}\n".format(i))
        optimized_geometry_file.write("{}\n".format(opt_molecule.to_xyz(fragments=[i])))

# STEP 1: calculate energies in the training set.
if args.calculate_training_set:
    calculator.fill_energies(settings_file_path,
                             args.training_set_input_path,
                             monomer_settings_paths,
                             optimized_geometry_paths,
                             args.training_set_output_path,
                             args.model_method,
                             args.model_basis,
                             args.model_counterpoise_correction)
else:
    args.training_set_output_path = args.training_set_input_path

# STEP 2: calculate charges, polarizabilities, and c6 constants!

if args.calculate_properties:
    potential_fitting.generate_fitting_config_file(settings_file_path, args.properties_path, geo_paths=optimized_geometry_paths, config_1b_paths=args.config_1b_paths, config_2b_paths=args.config_2b_paths)