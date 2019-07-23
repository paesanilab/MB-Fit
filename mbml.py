import sys, os, argparse
import potential_fitting
from potential_fitting import calculator
from potential_fitting.utils import SettingsReader, system
from potential_fitting.molecule import Molecule
from potential_fitting.exceptions import FunctionNotImplementedError, InvalidValueError

def string_to_boolean(string):
    if string.lower() == "true":
        return True
    elif string.lower() == "false":
        return False
    else:
        raise InvalidValueError("boolean", string, "'True' or 'False'")

# check arguments!

parser = argparse.ArgumentParser(description='Welcome to mbml! The python library for the generation of potential '
                                             'energy functions using permutationally invariant polynomials and QM methods.')

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
parser.add_argument('--properties_path', '-pp', dest='properties_path',
                                                  type=str, required=True,
                                                  help='File path to where the ".ini" file containing the properties'
                                                       'of the optimized geometry is to be written.')
parser.add_argument('--1b_config_paths', '-c1', nargs='+', dest='config_1b_paths', type=str, required=False,
                    help='List of file paths to ".ini" files containing 1b properties for each monomer. Only required '
                         'if your molecule has 2 or more fragments')
parser.add_argument('--2b_config_paths', '-c2', nargs='+', dest='config_2b_paths', type=str, required=False,
                    help='List of file paths to ".ini" files containing 2b properties for each dimer. Only required '
                         'if your molecule has 3 or more fragments')
parser.add_argument('--poly_directory', '-pd', dest='poly_directory_path', type=str, required=False, default=None,
                    help='File path to directory in which to store all polynomial files.')
parser.add_argument('--ttm_directory', '-td', dest='ttm_directory_path', type=str, required=False, default=None,
                    help='File path to directory in which to store all files for the ttm fit.')
parser.add_argument('--poly_fit_directory', '-fd', dest='poly_fit_directory_path', type=str, required=False, default=None,
                    help='File path to directory in which to store all files for the polynomial fit.')

parser.add_argument('--eval_script', '-es', dest='eval_script_path', type=str, required=False, default=None,
                    help='File path to file to write the script used to evaluate a configuration with the result of a fit.')
parser.add_argument('--fit_params', '-fp', dest='fit_params_path', type=str, required=False, default=None,
                    help='File path to file to write the paramaters used with the eval script to evaluate a configuration using the fit.')

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
parser.add_argument('--use_mbpol', '-um', dest='use_mbpol', type=str, required=True,
                    help='Comma delimited list of whether each fragment is water that should be modeled with lone pairs via mbpol.')

# The QM model

parser.add_argument('--method', '-mm', dest='model_method', type=str, required=False,
                    help='QM method to use for training set calculations.')
parser.add_argument('--basis', '-mb', dest='model_basis', type=str, required=False,
                    help='QM basis set to use for training set calculations.')
parser.add_argument('--counterpoise', '-mc', dest='model_counterpoise_correction', type=string_to_boolean, required=False,
                    help='If True, the counterpoise correction will be used for training set calculations.')

parser.add_argument('--properties_method', '-pm', dest='properties_method', type=str, required=False, default = 'wb97m-v',
                    help='QM method to use for calculating optimized geoemtry properties.')
parser.add_argument('--properties_basis', '-pb', dest='properties_basis', type=str, required=False, default = 'aug-cc-pvtz',
                    help='QM basis set to use for calculating optimized geometry properties.')

# Library options

parser.add_argument('--code', '-c', dest='code', type=str, required=True,
                    help='What code to use for calculations, options are "psi4" or "qchem".')
parser.add_argument('--num_threads', '-nt', dest='num_threads', type=str, required=False, default=1,
                    help='Number of threads to use for operations that support multithreading.')

# Fitting options

parser.add_argument('--poly_order', '-po', dest='poly_order', type=int, required=False,
                    help='Degree of polynomial to generate / use.')

skip_ttm_fit = parser.add_mutually_exclusive_group(required=False)
skip_ttm_fit.add_argument('--num_ttm_fits', '-ntf', dest='num_ttm_fits', type=int, required=False,
                    help='perform this many ttm fits and choose the best one.')
skip_poly_fit = parser.add_mutually_exclusive_group(required=True)
skip_poly_fit.add_argument('--num_poly_fits', '-npf', dest='num_poly_fits', type=int, required=False,
                    help='perform this many polynomail fits and choose the best one.')

# Optional arguments to skip various parts of the procedure.

skip_training_set_calculations_group.add_argument('--skip_training_set_calculations', '-stc', dest='calculate_training_set', required=False, action='store_false', default=True,
                    help='If included, then the input training set will be assumed to include calculated energies already, so no'
                         'calculations will be performed to put energies into the training set.')
parser.add_argument('--skip_properties_calculations', '-spc', dest='calculate_properties', required=False, action='store_false', default=True,
                    help='If included, then no calculation will be run to find the optimized geometry\'s properties. Instead, the properties path'
                         'specified with "--properties_path" will be assumed to already contain the properties.')
parser.add_argument('--skip_polynomial_generation', '-spg', dest='generate_polynomials', required=False, action='store_false', default=True,
                    help='If included, then no polynomials will be generated. Instead, the polynomials directory'
                         'specified with "--poly_directory" will be assumed to already contain all required polynomial files.')
skip_ttm_fit.add_argument('--skip_ttm_fit', '-stf', dest='perform_ttm_fit', required=False, action='store_false', default=True,
                    help='If included, then no ttm fit will be performed. Instead, the ttm directory'
                         'specified with "--ttm_directory" will be assumed to already contain all required files. Additionally, '
                         'if --skip_ttm_fit is specified, --skip_propertoes_calculations must also be specified.')
skip_poly_fit.add_argument('--skip_poly_fit', '-spf', dest='perform_poly_fit', required=False, action='store_false', default=True,
                    help='If included, then no poly fit will be performed. Instead, the poly fit directory'
                         'specified with "--poly_fit_directory" will be assumed to already contain all required files.')

args = parser.parse_args()

if not args.generate_polynomials and args.poly_directory_path is None:
    parser.error("Because --skip_polynomial_generation is specified, --poly_directory must be specified.")
if not args.perform_ttm_fit and args.ttm_directory_path is None:
    parser.error("Because --skip_ttm_fit is specified, --ttm_directory must be specified.")
if not args.perform_poly_fit and args.poly_fit_directory_path is None:
    parser.error("Because --skip_poly_fit is specified, --poly_fit_directory must be specified.")
if args.perform_poly_fit and args.eval_script_path is None:
    parser.error("Because --skip_poly_fit is not specified, --eval_script must be specified.")
if args.perform_poly_fit and args.fit_params_path is None:
    parser.error("Because --skip_poly_fit is not specified, --fit_params must be specified.")
if (args.generate_polynomials or args.perform_poly_fit) and args.poly_order is None:
    parser.error("Because either --skip_polynomial_generation or --skip_poly_fit are not specified, you must specify the order of the polynomial with --poly_order")
if args.calculate_training_set:
    if args.model_method is None:
        parser.error("Because --skip_training_set_calculations is not specified, you must specify the method to use for training set calculations with --method.")
    if args.model_basis is None:
        parser.error("Because --skip_training_set_calculations is not specified, you must specify the basis to use for training set calculations with --basis.")
    if args.number_of_fragments > 1 and args.model_counterpoise_correction is None:
        parser.error("Because number of fragments is greater than 1 and --skip_training_set_calculations is not specified, you must specify the counterpoise corretion to use for training set calculations with --counterpoise.")

if args.number_of_fragments == 1 and args.model_counterpoise_correction:
    parser.error("Because number of fragments is 1, counterpoise correction cannot be True")

if not args.perform_ttm_fit and args.calculate_properties:
    parser.error("Because --skip_ttm_fit is specified, --skip_properties_calculations must also be specified.")
if args.number_of_fragments > 1 and args.perform_ttm_fit and args.num_ttm_fits is None:
    parser.error("Because number of fragments is greater than 1 and --skip_ttm_fit is not specified, you must specify the number of ttm fits to perform with --num_ttm_fits")

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
settings.set('molecule', 'use_mbpol', args.use_mbpol)

settings.set('model', 'method', args.model_method)
settings.set('model', 'basis', args.model_basis)
settings.set('model', 'cp', args.model_counterpoise_correction)

settings.set('energy_calculator', 'code', args.code)
settings.set('qchem', 'num_threads', args.num_threads)
settings.set('psi4', 'num_threads', args.num_threads)

settings.set('config', 'method', args.properties_method)
settings.set('config', 'basis', args.properties_basis)

molecule_in = "_".join(settings.get('molecule', 'symmetry').split(','))

temp_file_path = os.path.join(settings.get('files', 'log_path'), "temp_files")
settings_file_path = os.path.join(temp_file_path, "settings.ini")

if args.poly_directory_path is None:
    args.poly_directory_path = os.path.join(temp_file_path, "poly_directory")
if args.ttm_directory_path is None:
    args.ttm_directory_path = os.path.join(temp_file_path, "ttm_directory")
if args.poly_fit_directory_path is None:
    args.poly_fit_directory_path = os.path.join(temp_file_path, "poly_fit_directory")

poly_in_path = os.path.join(args.poly_directory_path, "poly.in")

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

optimized_geometry_paths = [os.path.join(temp_file_path, "opt_geo{}.xyz".format(i)) for i in range(args.number_of_fragments)]

opt_molecule = Molecule.read_xyz_path_direct(args.optimized_geometry_path, settings)

for i, optimized_geometry_path in enumerate(optimized_geometry_paths):
    with open(optimized_geometry_path, "w") as optimized_geometry_file:
        optimized_geometry_file.write("{}\n".format(opt_molecule.get_fragments()[i].get_num_atoms()))
        optimized_geometry_file.write("optimized geometry for fragment {}\n".format(i))
        optimized_geometry_file.write("{}\n".format(opt_molecule.to_xyz(fragments=[i])))

# STEP 1: calculate energies in the training set.
if args.calculate_training_set:
    system.format_print("Finding energies of training set...", bold=True, italics=True, color=system.Color.YELLOW)
    calculator.fill_energies(settings_file_path,
                             args.training_set_input_path,
                             monomer_settings_paths,
                             optimized_geometry_paths,
                             args.training_set_output_path,
                             args.model_method,
                             args.model_basis,
                             args.model_counterpoise_correction)
    system.format_print("Training set energies calculated successfully!", bold=True, italics=True, color=system.Color.GREEN)
else:
    system.format_print("Training set energies already calculated, no need to caclulate them.", bold=True, italics=True, color=system.Color.BLUE)
    args.training_set_output_path = args.training_set_input_path

# STEP 2: calculate charges, polarizabilities, and c6 constants!

if args.calculate_properties:
    system.format_print("Finding properties of optimized geometry...", bold=True, italics=True, color=system.Color.YELLOW)
    potential_fitting.generate_fitting_config_file(settings_file_path, args.properties_path, geo_paths=optimized_geometry_paths, config_1b_paths=args.config_1b_paths, config_2b_paths=args.config_2b_paths)
    system.format_print("Optimized properties calculated successfully!", bold=True, italics=True, color=system.Color.GREEN)
else:
    system.format_print("Optimized properties already calculated, no need to caclulate them.", bold=True, italics=True, color=system.Color.BLUE)

# STEP 3 generate polynomials

if args.generate_polynomials:
    system.format_print("Generating polynomials...", bold=True, italics=True, color=system.Color.YELLOW)
    potential_fitting.generate_poly_input(settings_file_path, molecule_in, poly_in_path)
    potential_fitting.generate_polynomials(settings_file_path, poly_in_path, args.poly_order, args.poly_directory_path)
    potential_fitting.execute_maple(settings_file_path, args.poly_directory_path)
    system.format_print("Polynomial generation successful!", bold=True, italics=True, color=system.Color.GREEN)
else:
    system.format_print("Polynomials already generated, no need to generate them.", bold=True, italics=True, color=system.Color.BLUE)

# STEP 4: for 2B+ run the TTM fit!

if opt_molecule.get_num_fragments() == 1:
    system.format_print("Molecule only has one fragment, no need to perform TTM fit.", bold=True, italics=True, color=system.Color.BLUE)
elif args.perform_ttm_fit:
    system.format_print("Performing TTM fit...", bold=True, italics=True, color=system.Color.YELLOW)
    if opt_molecule.get_num_fragments() == 2:
        potential_fitting.generate_2b_ttm_fit_code(settings_file_path, args.properties_path, molecule_in, args.ttm_directory_path)
        potential_fitting.compile_fit_code(settings_file_path, args.ttm_directory_path)
        potential_fitting.fit_2b_ttm_training_set(settings_file_path, os.path.join(args.ttm_directory_path, "fit-2b-ttm"), args.training_set_output_path, args.ttm_directory_path, args.properties_path, args.num_ttm_fits)
    else:
        raise FunctionNotImplementedError("ttm fits for 3b+")
    system.format_print("TTM fit successful!", bold=True, italics=True, color=system.Color.GREEN)
else:
    system.format_print("TTM fit already perfomed, no need to run it again.", bold=True, italics=True, color=system.Color.BLUE)

# STEP 5: run the polynomial fit!

if args.perform_poly_fit:
    system.format_print("Performing Polynomial fit...", bold=True, italics=True, color=system.Color.YELLOW)
    if opt_molecule.get_num_fragments() == 1:
        potential_fitting.generate_1b_fit_code(settings_file_path, args.properties_path, molecule_in, poly_in_path, args.poly_directory_path, args.poly_order, args.poly_fit_directory_path)
        potential_fitting.compile_fit_code(settings_file_path, args.poly_fit_directory_path)
        system.call("mv", os.path.join(args.poly_fit_directory_path, "eval-1b"), args.eval_script_path)
        potential_fitting.fit_1b_training_set(settings_file_path, os.path.join(args.poly_fit_directory_path, "fit-1b"), args.training_set_output_path, args.poly_fit_directory_path, args.fit_params_path, args.num_poly_fits)
    elif opt_molecule.get_num_fragments() == 2:
        potential_fitting.generate_2b_fit_code(settings_file_path, args.properties_path, poly_in_path, args.poly_directory_path, args.poly_order, args.poly_fit_directory_path)
        potential_fitting.compile_fit_code(settings_file_path, args.poly_fit_directory_path)
        system.call("mv", os.path.join(args.poly_fit_directory_path, "eval-2b"), args.eval_script_path)
        potential_fitting.fit_2b_training_set(settings_file_path, os.path.join(args.poly_fit_directory_path, "fit-2b"), args.training_set_output_path, args.poly_fit_directory_path, args.fit_params_path, args.num_poly_fits)
    else:
        raise FunctionNotImplementedError("polynomial fits for 3b+")
    system.format_print("Polynomial fit successful!", bold=True, italics=True, color=system.Color.GREEN)
else:
    system.format_print("Polynomial fit already perfomed, no need to run it again.", bold=True, italics=True, color=system.Color.BLUE)


# STEP 6: test the fit!
