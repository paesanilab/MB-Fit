# absolute module imports
from potential_fitting import calculator
from potential_fitting.molecule import xyz_to_molecules
from potential_fitting.utils import SettingsReader, files, system
from potential_fitting.calculator import Model

def generate_normal_modes(settings_path, opt_path, normal_modes_path, method, basis, arguments={}):
    """
    Generates the a normal modes input file from an optimized geometry.

    Args:
        settings_path       - Local path to the ".ini" file containing all relevant settings.
        opt_path            - Local path to the ".xyz" file to read the optimized geometry from.
        normal_modes_path   - Local path to the ".dat" file to write the normal modes to.
        arguments           - Dictionary of extra arguments to be passed to the QM code doing the calculation.

    Returns:
        The null dimension of the input molecule.
    """
    system.format_print("Beginning normal modes calculation of {}...".format(opt_path),
            bold=True, color=system.Color.YELLOW)
    
    settings = SettingsReader(settings_path)

    calc = calculator.get_calculator(settings_path)
    model = Model(method, basis)

    # parse the optimized geometry
    molecule = xyz_to_molecules(opt_path, settings)[0]
    
    # calculate the normal modes
    normal_modes, frequencies, red_masses, log_path = calc.find_frequencies(molecule, model, arguments=arguments)

    # calculate dim null
    dim_null = 3 * molecule.get_num_atoms() - len(normal_modes)
    
    # write the normal modes to the output file
    write_normal_modes(settings_path, normal_modes, frequencies, red_masses, normal_modes_path)

    system.format_print("Normal modes calculation complete! Normal modes in {}.".format(normal_modes_path),
            bold=True, color=system.Color.GREEN)

def write_normal_modes(settings_path, normal_modes, frequencies, red_masses, normal_modes_path):
    """
    Writes the info returned by calculator.frequencies() to the file normal_modes_path.

    Args:
        settings_path       - Local path to the ".ini" file containing all relevant settings.
        normal_modes        - The normal modes to write.
        frequencies         - The frequencies to write.
        red_masses          - The reduced_masses to write.
        normal_modes_path   - Local path to the file to write the normal modes to.
    """

    settings = SettingsReader(settings_path)

    # formatter strings
    norm_formatter = "{:> 12.4f}"
    freq_formatter = "{:.2f}"
    mass_formatter = "{:.6f}"

    normal_out = ""
    
    num_modes = len(frequencies)

    # for each normal mode
    for i in range(1, 1 + num_modes):
        index = i - 1
        
        # get the geometry of this normal mode
        geo = normal_modes[index]
        
        normal_out += "normal mode: " + str(i) + "\n"
        normal_out += "frequency = " + freq_formatter.format(frequencies[index]) + "\n"
        normal_out += "reduced mass = " + mass_formatter.format(red_masses[index]) + "\n"
        
        # for each atom in the molecule
        for atom in range(len(geo)):
            # write this atom's x, y, and z for this normal mode, setting any values below .000001 to 0.
            normal_out += norm_formatter.format(0.0 if abs(float(geo[atom][0])) < 1e-6 else float(geo[atom][0])) + "\t"
            normal_out += norm_formatter.format(0.0 if abs(float(geo[atom][1])) < 1e-6 else float(geo[atom][1])) + "\t"
            normal_out += norm_formatter.format(0.0 if abs(float(geo[atom][2])) < 1e-6 else float(geo[atom][2])) + "\n"
        
        normal_out += "\n"

    normal_modes_path = files.init_file(normal_modes_path, files.OverwriteMethod.get_from_settings(settings))

    # write the normal modes to the output file
    with open(normal_modes_path, 'w') as norm_file:
        norm_file.write(normal_out)
