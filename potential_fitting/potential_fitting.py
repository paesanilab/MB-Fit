# external package imports
import os, subprocess, contextlib, random

# local module imports
from .utils import SettingsReader, files, system
from . import configurations, database, polynomials, fitting
from .database import Database

def optimize_geometry(settings_path, unopt_geo_path, opt_geo_path):
    """
    Optimizes the geometry of the given molecule.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        unopt_geo_path      - Local path to the file to read the unoptimized geoemtry from.
        opt_geo_path        - Local path to the file to write the optimized geometry to.

    Returns:
        None.
    """

    configurations.optimize_geometry(settings_path, unopt_geo_path, opt_geo_path)

def generate_normal_modes(settings_path, opt_geo_path, normal_modes_path):
    """
    Generates the normal modes for the given molecule.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        opt_geo_path        - Local path to the file to read the optimized geometry from.
        normal_modes_path   - Local path to the file to write the normal modes to.

    Returns:
        Null dimension of normal modes.
    """
    
    dim_null = configurations.generate_normal_modes(settings_path, opt_geo_path, normal_modes_path)

    return dim_null

def generate_1b_configurations(settings_path, opt_geo_path, normal_modes_path, configurations_path,
        number_of_configs = 100, seed = None):
    """
    Generates 1b configurations for a given monomer from a set of normal modes.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        opt_geo_path        - Local path to the file to read optimized geometry from.
        normal_modes_path   - Local path to the file to read normal modes from.
        dim_null            - The null dimension of this molecule, see generate_normal_modes().
        config_path         - Local path to the file to write configurations to.
        number_of_configs   - Number of configurations to generate
        seed                - The same seed with the same molecule and normal modes will always generate the same
                configurations.

    Returns:
        None.
    """

    if seed is None:
        seed = random.randint(-10000000, 10000000)

    configurations.generate_1b_configurations(settings_path, opt_geo_path, normal_modes_path, configurations_path,
            number_of_configs, seed = seed)

def generate_2b_configurations(settings_path, geo1_path, geo2_path, number_of_configs, configurations_path, 
        min_distance = 1, max_distance = 5, min_inter_distance = 0.8, progression = False, use_grid = False, 
        step_size = 0.5, seed = None):
    """
    Generates 2b configurations for a given dimer.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        geo1_path           - Local path to read the first optimized geometry from.
        geo2_path           - Local path to read the second optimized geometry from.
        number_of_configs   - The target number of configurations to generate. If max_distance is set too low or 
                min_inter_distance is set too high, then less configurations may be generated.
        configurations_path - Local path to the file in to write the configurations to.
        min_distance        - The minimum distance between the centers of mass of the two molecules in any
                configuration.
        max_distance        - The maximum distance between the centers of mass of the two molecules in any
                configuration.
        min_inter_distance  - Minimum intermolecular distance in any config is the sum of two atoms vdw radii * this
                value.
        progression         - If True, use a linear progression for distance between atoms, otherwise random
                progression.r
        use_grid            - If False, configurations are space roughly evenly between min_distance and max_distance.
                If True, then configurations are placed at intervals along this distance based on step_size.
        step_size           - If use_grid is True, then this dictates the distance of the spacing interval used to
                place the centers of masses of the molecules. Otherwise, this parameter has no effect.
        seed                - The same seed will generate the same configurations.

    Returns:
        None.
    """

    if seed is None:
        seed = random.randint(-10000000, 10000000)

    configurations.generate_2b_configurations(geo1_path, geo2_path, number_of_configs, configurations_path,
            min_distance, max_distance, min_inter_distance, progression, use_grid, step_size, seed)

def init_database(settings_path, database_path, configurations_path, *tags):
    """
    Creates a database from the given configuration .xyz files. Can be called on a new database
    to create a new database, or an existing database to add more energies to be calculated

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        database_path       - Local path to the database file to add the configurations to. ".db" will automatically
                be added to the end if it does not already end in ".db".
        configurations_path - Local path to a single .xyz file or a directory containing ".xyz" files. If this argument
                is a directory, any non .xyz files will be ignored.
        tags                - Mark the new configurations with these tags.

    Returns:
        None.
    """

    database.initialize_database(settings_path, database_path, configurations_path, *tags)

def fill_database(settings_path, database_path):
    """
    Goes through all the uncalculated energies in a database and calculates them. Will take a while. May be interrupted
    and restarted.
    
    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        database_path       - Local path to the database file containing uncaculated energies. ".db" will
                automatically be added to the end if it does not already end in ".db".

    Returns:
        None.
    """

    database.fill_database(settings_path, database_path)

def generate_1b_training_set(settings_path, database_path, training_set_path, molecule_name, *tags, method = "%", basis = "%",
            cp = "%", e_min = 0, e_max = float('inf')):
    """
    Generates a 1b training set from the energies inside a database.

    Specific method, bnasis, and cp may be specified to only use energies calculated
    with a specific model.

    '%' can be used to stand in as a wild card, meaning any method/basis/cp will be used in the training set.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        database_path       - Local path to the database file containing energies to put in the training set. ".db"
                will automatically be added to the end if it does not already end in ".db".
        training_set_path   - Local path to the file to write the training set to.
        molecule_name       - The name of the moelcule to generate a training set for.
        tags                - Only use energies marked with one or more of these tags.
        method              - Only use energies calcualted by this method.
        basis               - Only use energies calculated in this basis.
        cp                  - Only use energies calculated with this cp. Note that counterpoise correct has no
                effect on 1b energies.
        e_min               - Minimum (inclusive) energy of any config to include in this training set.
        e_max               - Maximum (exclusive) energy of any config to include in this training set.

    Returns:
        None.
    """

    database.generate_1b_training_set(settings_path, database_path, training_set_path, molecule_name,
            method, basis, cp, *tags, e_min = e_min, e_max = e_max)

def generate_2b_training_set(settings_path, database_path, training_set_path, monomer1_name, monomer2_name, *tags,
        method = "%", basis = "%", cp = "%", e_bind_max = float('inf'), e_mon_max = float('inf')):
    """
    Generates a 2b training set from the energies inside a database.

    Specific method, basis, and cp may be specified to only use energies calculated
    with a specific model.

    '%' can be used to stand in as a wild card, meaning any method/basis/cp will be used in the training set.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        database_path       - Local path to the database file containing energies to put in the training set. ".db"
                will automatically be added to the end if it does not already end in ".db".
        training_set_path   - Local path to the file to write the training set to.
        monomer1_name       - The Name of first monomer in the dimer.
        monomer2_name       - The Name of the second monomer in the dimer.
        tags                - Only use energies marked with one or more of these tags.
        method              - Only use energies calcualted by this method.
        basis               - Only use energies calculated in this basis.
        cp                  - Only use energies calculated with this cp.
        e_bind_max          - Maximum binding energy allowed
        e_mon_max           - Maximum monomer deformation energy allowed

    Returns:
        None.
    """
    
    database.generate_2b_training_set(settings_path, database_path, training_set_path, monomer1_name, monomer2_name,
            method, basis, cp, *tags, e_bind_max = e_bind_max, e_mon_max = e_mon_max)

def generate_poly_input(settings_path, molecule_in, in_file_path):
    """
    Generates an input file for polynomial generation.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        molecule_in         - String idicating symmetry of molecule, ie "A1B2_A1B2" for (CO2)2
        in_file_path        - Local path to the file to write the polynomial input to.

    Returns:
        None.
    """

    polynomials.generate_input_poly(settings_path, molecule_in, in_file_path)

def generate_poly_input_from_database(settings_path, database_path, molecule_name, in_file_path):
    """
    Generates an input file for polynomial generation.
    Looks in a database to find the symmetry and creates a file in the given directory.

    If the symmetry is A1B2, then the file A1B2.in containing polynomial generation input will be created inside
    the poly_directory_path directory.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        database_path       - Local path to the database file containing the molecular symmetry. ".db" will
                automatically be added to the end if it does not already end in ".db".
        molecule_name       - The name of the molecule to generate a polynomial generation input file for. At least one
                instance of this molecule must be in the database.
        in_file_path        - Local path to the file to write the polynomial input to.

    Returns:
        None.
    """

    with Database(database_path) as database:
        symmetry = database.get_symmetry(molecule_name)

        generate_poly_input(settings_path, symmetry, in_file_path)

def generate_polynomials(settings_path, poly_in_path, order, poly_dir_path):
    """
    Generates polynomial input for maple and some ".cpp" and ".h" polynomial files.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        poly_in_path        - Local path to the file to read the polynomial input from. Name of file should be in
                the format "A1B2.in". It is ok to have extra directories prior to file name. For example:
                "thisplace/thatplace/A3.in".
        order               - The order of the polynomial to generate.
        poly_dir_path       - Local path to the directory to write the polynomial files in.

    Returns:
        None.
    """

    polynomials.generate_poly(settings_path, poly_in_path, order, poly_dir_path)

def execute_maple(settings_path, poly_dir_path):
    """
    Runs maple on the polynomial files in the specified directory to turn them into actual cpp files.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        poly_directory      - Local path to the directory to read the ".maple" files from and write the ".cpp" files
                to.

    Returns:
        None.
    """

    this_file_path = os.path.dirname(os.path.abspath(__file__))

    original_dir = os.getcwd()

    os.chdir(poly_dir_path)

    # clear any old files, because for some reason maple appends to existing files instead of clobbering
    with contextlib.suppress(FileNotFoundError):
        os.remove("poly-grd.c")
        os.remove("poly-nogrd.c")
        os.remove("poly-grd.cpp")
        os.remove("poly-nogrd.cpp")
    
    system.call("maple", "poly-grd.maple")
    system.call("maple", "poly-nogrd.maple")

    with open("poly-grd.c", "r") as in_file, open("poly-grd.cpp", "w") as out_file:
        system.call("clean-maple-c.pl", in_file = in_file, out_file = out_file)
    with open("poly-nogrd.c", "r") as in_file, open("poly-nogrd.cpp", "w") as out_file:
        system.call("clean-maple-c.pl", in_file = in_file, out_file = out_file)

    os.chdir(original_dir)

def generate_fit_config(settings_path, molecule_in, config_path, *opt_geometry_paths, distance_between = 20):
    """
    Generates the config file needed to perform a fit from the optimized geometries of up to 3 monomers.

    Qchem is required for this step to work.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        molecule_in         - A String of fromat "A1B2". Same as poly_in_path but without ".in".
        config_path         - Local path to file to write the config file to.
        opt_geometry_paths  - Local paths to the optimized geometries to include in this fit config, should be 1 to 3
                (inclusive) of them.
        distance_between    - The Distance between each geometry in the qchem calculation. If the qchem calculation
                does not converge, try different values of this.

    Returns:
        None.
    """

    fitting.make_config(settings_path, molecule_in, config_path, *opt_geometry_paths,
            distance_between = distance_between)
    
def generate_1b_fit_code(settings_path, config_path, molecule_in, poly_in_path, poly_dir_path, order, fit_dir_path):
    """
    Generates the fit code based on the polynomials and config file for a monomer.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        config_path         - Local path to the monomer config file to read config info from.
        poly_in_path        - Local path to the file to read the polynomial input from. Name of file should be in
                the format A1B2.in, it is ok to have extra directories prior to file name (thisplace/thatplace/A3.in).
        poly_dir_path       - Local path to the directory where the polynomial ".h" and ".cpp" files are files are.
        order               - The order of the polynomial to in poly_dir_path.
        fit_dir_path        - Local path to the directory to write the fit code in.

    Returns:
        None.
    """

    fitting.prepare_1b_fitting_code(config_path, molecule_in, poly_in_path, poly_dir_path, order, fit_dir_path)

def generate_2b_ttm_fit_code(settings_path, config_path, molecule_in, fit_dir_path):
    """
    Generates the fit TTM fit code for a dimer.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        config_path         - Local path to the dimer config file to read config info from.
        molecule_in         - A String of fromat "A1B2". Same as poly_in_path but without ".in".
        fit_dir_path        - Local path to the directory to write the ttm fit code in.

    Returns:
        None.
    """

    this_file_path = os.path.dirname(os.path.abspath(__file__))

    original_dir = os.getcwd()

    files.init_directory(fit_dir_path)

    os.chdir(fit_dir_path)

    codes_path = os.path.join(this_file_path, "..", "codes", "2b-codes")
    
    # FOR SOME REASON THE SECOND LINE WORKS BUT NOT THE FIRST, WHY? WHO KNOWS.
    # system.call("cp", os.path.join(codes_path, "template/*"), ".")
    os.system("cp " + os.path.join(codes_path, "template", "*") + " .")

    ttm_script_path = os.path.join(codes_path, "get_2b_TTM_codes.py")

    system.call("python", ttm_script_path, "{}/{}".format(original_dir, settings_path), "{}/{}".format(original_dir, config_path), molecule_in)

    os.chdir(original_dir)   
 
def generate_2b_fit_code(settings_path, config_path, poly_in_path, poly_path, poly_order, fit_dir_path):
    """
    Generates the fit code based on the polynomials for a monomer

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        config_path         - Local path to the dimer config file.
        poly_in_path        - Local path to the the A3B2.in type file to read polynomial input from.
        poly_path           - Local path to directory where polynomial files are.
        poly_order          - The order of the polynomial in poly_path.
        fit_dir_path        - Local path to directory to generate fit code in.

    Returns:
        None.
    """

    files.init_directory(fit_dir_path)

    if not os.path.isdir(fit_dir_path):
        os.mkdir(fit_dir_path)

    fitting.prepare_2b_fitting_code(settings_path, config_path, poly_in_path, poly_path, poly_order, fit_dir_path)
 
def compile_fit_code(settings_path, fit_dir_path):
    """
    Compiles the fit code in the given directory.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fit_dir_path        - Local path to the directory to read uncompiled fit code from and write compiled fit code
                to.

    Returns:
        None.
    """

    original_dir = os.getcwd()

    os.chdir(fit_dir_path)

    system.call("make", "clean")
    system.call("make")

    os.chdir(original_dir)

def fit_1b_training_set(settings_path, fit_code_path, training_set_path, fit_dir_path, fitted_nc_path, num_fits = 10):
    """
    Fits a given 1b training set using a given 1b fit code.

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fit_code_path       - Local path to the fit code with which to fit the training set. Generated in fit_dir_path
                by compile_fit_code.
        training_set_path   - Local path to the file to read the training set from.
        fit_dir_path        - Local path to the directory to write the ".cdl" and ".dat" files created by this fit.
        fitted_nc_path      - Local path to the file to write the final fitted ".nc" to.
        num_fits            - Performs this many fits and then gives only the best one.

    Returns:
        None
    """

    settings = SettingsReader(settings_path)

    # init the required files
    files.init_directory(fit_dir_path)
    best_fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "mb", "best-fit.log"))
    fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "mb" "fit.log"))

    # keeps tracks of the number of attempts to generate a fit
    attempts = 1

    # generate an initial fit
    with open(best_fit_log_path, "w") as best_fit_log:
        system.call(fit_code_path, training_set_path, out_file = best_fit_log)
        os.rename("fit-1b.cdl", "best-fit-1b.cdl")
        os.rename("fit-1b-initial.cdl", "best-fit-1b-initial.cdl")
        os.rename("correlation.dat", "best-correlation.dat")

    with open(best_fit_log_path, "r") as best_fit_log:
        best_log_lines = best_fit_log.readlines()

    best_rmsd = float(best_log_lines[-6].split()[2])

    print("Completed first fit with rmsd {}.\n".format(best_rmsd))

    while(attempts < num_fits):

        # generate a new fit
        with open(fit_log_path, "w") as fit_log:
            system.call(fit_code_path, training_set_path, out_file = fit_log)
  
        with open(fit_log_path, "r") as fit_log:
            log_lines = fit_log.readlines()

        rmsd = float(log_lines[-6].split()[2])

        print("Completed fit number {} with rmsd {}.".format(attempts, rmsd))

        print("Current best fit has rmsd {}.".format(best_rmsd))

        # if the new fit is better than the old fit, replace the best log and best cdl files
        if rmsd < best_rmsd:

            print("Replaced previous best fit with most recent one.")

            os.rename(fit_log_path, best_fit_log_path)
            os.rename("fit-1b.cdl", "best-fit-1b.cdl")
            os.rename("fit-1b-initial.cdl", "best-fit-1b-initial.cdl")
            os.rename("correlation.dat", "best-correlation.dat")

            best_rmsd = rmsd
            
        attempts += 1

        print("\n")

    # remove the most recent fit file
    try:
        os.remove(fit_log_path)
        os.remove("fit-1b.cdl")
        os.remove("fit-1b-initial.cdl")
        os.remove("correlation.dat")
    # in the case that there is no most recent fit file because the last fit was the best fit, do nothing
    except FileNotFoundError:
        pass

    system.call("ncgen", "-o", fitted_nc_path, "best-fit-1b.cdl")

    os.rename("best-fit-1b.cdl", os.path.join(fit_dir_path, "fit-1b.cdl"))
    os.rename("best-fit-1b-initial.cdl", os.path.join(fit_dir_path, "fit-1b-initial.cdl"))
    os.rename("best-correlation.dat", os.path.join(fit_dir_path, "correlation.dat"))

def fit_2b_ttm_training_set(settings_path, fit_code_path, training_set_path, fit_dir_path, config_path, num_fits = 10):
    """
    Fits a given 2b training set using a given 2b ttm fit code

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fit_code_path       - Local path to the fit code with which to fit the training set. Generated in fit_dir_path
                by compile_fit_code.
        training_set_path   - Local path to the file to read the training set from.
        fit_dir_path        - the directory where the fit log and other files created by the fit go
        config_path         - Local path to the ".ini" file to write A6 and d6 constants to.
        num_fits            - Performs this many fits and then gives only the best one.

    Returns:
        None
    """

    settings = SettingsReader(settings_path)

    files.init_directory(fit_dir_path)

    best_fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "ttm", "best_fit.log"))
    fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "ttm", "fit.log"))

    attempts = 1;
    with open(best_fit_log_path, "w") as best_fit_log:
        system.call(fit_code_path, training_set_path, out_file = best_fit_log)
        os.rename("individual_terms.dat", "best-individual_terms.dat")
        os.rename("ttm-params.txt", "best-ttm-params.txt")
        os.rename("correlation.dat", "best-correlation.dat")

    with open(best_fit_log_path, "r") as best_fit_log:
        best_log_lines = best_fit_log.readlines()

    best_rmsd = float(best_log_lines[-4].split()[2])

    print("Completed first fit with rmsd {}.\n".format(best_rmsd))

    while(attempts < num_fits):

        with open(fit_log_path, "w") as fit_log:
            system.call(fit_code_path, training_set_path, out_file = fit_log)
        
        with open(fit_log_path, "r") as fit_log:
            log_lines = fit_log.readlines()

        rmsd = float(log_lines[-4].split()[2])

        print("Completed fit number {} with rmsd {}.".format(attempts, rmsd))

        print("Current best fit has rmsd {}.".format(best_rmsd))

        if rmsd < best_rmsd:

            print("Replaced previous best fit with most recent one.")

            os.rename(fit_log_path, best_fit_log_path)
            os.rename("individual_terms.dat", "best-individual_terms.dat")
            os.rename("ttm-params.txt", "best-ttm-params.txt")
            os.rename("correlation.dat", "best-correlation.dat")

            best_rmsd = rmsd
            

        attempts += 1

        print("\n")

    # remove the most recent fit file
    try:
        os.remove(fit_log_path)
        os.remove("individual_terms.dat")
        os.remove("ttm-params.txt")
        os.remove("correlation.dat")
    # in the case that there is no most recent fit file because the last fit was the best fit, do nothing
    except FileNotFoundError:
        pass

    os.rename("best-individual_terms.dat", os.path.join(fit_dir_path, "individual_terms.dat"))
    os.rename("best-ttm-params.txt", os.path.join(fit_dir_path, "ttm-params.txt"))
    os.rename("best-correlation.dat", os.path.join(fit_dir_path, "correlation.dat"))

    # read d6 and A constants from ttm output file
    with open(os.path.join(fit_dir_path, "ttm-params.txt"), "r") as ttm_file:
        A = [float(a) for a in ttm_file.readline().split()]
        d6 = [float(d) for d in ttm_file.readline().split()]

    # write d6 and A to the config.ini file
    lines = []
    with open(config_path, "r") as config_file:
        for line in config_file:
            lines.append(line)

    found_A = False
    with open(config_path, "w") as config_file:
        for line in lines:
            if line.startswith("A = "):
                found_A = True
            if not line.startswith("d6 = "):
                config_file.write(line)

            else:
                config_file.write("d6 = {}\n".format([[], [], d6]))

        if not found_A:
            config_file.write("A = {}\n".format([[], [], A]))

def fit_2b_training_set(settings_path, fit_code_path, training_set_path, fit_dir_path, fitted_nc_path, num_fits = 10):
    """
    Fits the fit code to a given training set

    Args:
        settings_path       - Local path to the file containing all relevent settings information.
        fit_code            - Local path to the code to fit.
        training_set        - Local path to the training set to fit the code to.
        fit_dir_path        - Local path to the directory where the .cdl and .dat files will end up.
        fitted_nc_path      - Local path to file to write final fitted ".nc" file to.
        num_fits            - Performs this many fits and then gives only the best one.

    Returns:
        None
    """

    settings = SettingsReader(settings_path)

    # init the required files
    files.init_directory(fit_dir_path)
    best_fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "mb", "best-fit.log"))
    fit_log_path = files.init_file(os.path.join(settings.get("files", "log_path"), "mb" "fit.log"))

    # keeps tracks of the number of attempts to generate a fit
    attempts = 1

    # generate an initial fit
    with open(best_fit_log_path, "w") as best_fit_log:
        system.call(fit_code_path, training_set_path, out_file = best_fit_log)
        os.rename("fit-2b.cdl", "best-fit-2b.cdl")
        os.rename("fit-2b-initial.cdl", "best-fit-2b-initial.cdl")
        os.rename("correlation.dat", "best-correlation.dat")

    with open(best_fit_log_path, "r") as best_fit_log:
        best_log_lines = best_fit_log.readlines()

    best_rmsd = float(best_log_lines[-6].split()[2])

    print("Completed first fit with rmsd {}.\n".format(best_rmsd))

    while(attempts < num_fits):

        # generate a new fit
        with open(fit_log_path, "w") as fit_log:
            system.call(fit_code_path, training_set_path, out_file = fit_log)
  
        with open(fit_log_path, "r") as fit_log:
            log_lines = fit_log.readlines()

        rmsd = float(log_lines[-6].split()[2])

        print("Completed fit number {} with rmsd {}.".format(attempts, rmsd))

        print("Current best fit has rmsd {}.".format(best_rmsd))

        # if the new fit is better than the old fit, replace the best log and best cdl files
        if rmsd < best_rmsd:

            print("Replaced previous best fit with most recent one.")

            os.rename(fit_log_path, best_fit_log_path)
            os.rename("fit-2b.cdl", "best-fit-2b.cdl")
            os.rename("fit-2b-initial.cdl", "best-fit-2b-initial.cdl")
            os.rename("correlation.dat", "best-correlation.dat")

            best_rmsd = rmsd
            
        attempts += 1

        print("\n")

    # remove the most recent fit file
    try:
        os.remove(fit_log_path)
        os.remove("fit-2b.cdl")
        os.remove("fit-2b-initial.cdl")
        os.remove("correlation.dat")
    # in the case that there is no most recent fit file because the last fit was the best fit, do nothing
    except FileNotFoundError:
        pass

    system.call("ncgen", "-o", fitted_nc_path, "best-fit-2b.cdl")

    os.rename("best-fit-2b.cdl", os.path.join(fit_dir_path, "fit-2b.cdl"))
    os.rename("best-fit-2b-initial.cdl", os.path.join(fit_dir_path, "fit-2b-initial.cdl"))
    os.rename("best-correlation.dat", os.path.join(fit_dir_path, "correlation.dat"))
