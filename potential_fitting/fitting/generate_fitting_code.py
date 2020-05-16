import sys, os, json
from potential_fitting.utils import SettingsReader
from potential_fitting.utils import constants, system
from potential_fitting.polynomials import MoleculeSymmetryParser
from . import utils_nb_fitting 
from . import file_writer_nb_fitting 


def generate_fitting_code(settings_path, config_path, degree, poly_in_path, poly_directory_path, use_direct):

    directcpp = poly_directory_path + "/poly-direct.cpp"
    directgradcpp = poly_directory_path + "/poly-grd-direct.cpp"
    name = poly_in_path
    poly_directory = poly_directory_path
    
    # In[ ]:
    
    settings = SettingsReader(settings_path)
    config = SettingsReader(config_path)
    
    
    ################################################################################
    ## MONOMER PROPERTIES ##########################################################
    ################################################################################
    
    # Define list of variables that are fictitious
    virtual_sites_poly = config.getlist("fitting", "virtual_site_labels")
    
    symmetry_parser = MoleculeSymmetryParser("_".join(settings.get("molecule", "symmetry").split(",")), virtual_sites=virtual_sites_poly)
    
    system.format_print("Generating fitcode for molecule with symmetry {}...".format(symmetry_parser.get_symmetry()),
                        bold=True, color=system.Color.YELLOW)
    
    # Obtain monomers from settings
    monomers_symmetries = list(symmetry_parser.get_fragment_symmetries())
    
    # Get number of atoms in each monomer.
    number_of_atoms = [int(i) for i in settings.get("molecule", "fragments").split(",")]
    number_of_monomers = len(number_of_atoms)
    
    # Set the number of electrostatic sites
    number_of_sites = [i for i in number_of_atoms]
    
    # Get which monomer should be mb-pol monomer (if any) 1 means yes, 0 means no.
    use_mbpol = [int(i) for i in settings.get("molecule", "use_mbpol").split(",")]
    
    system.format_print("Using mbpol for {} fragments.".format(sum(use_mbpol)),
                        italics=True)
    
    # Define if lone pairs are used based on monomer names
    # Update number of sites if needed
    use_lonepairs = [0]*number_of_monomers
    for i in range(number_of_monomers):
        if any(x in monomers_symmetries[i] for x in virtual_sites_poly):
            use_lonepairs[i] = 1
        if use_mbpol[i] != 0:
            number_of_sites[i] += 1
    
    system.format_print("{} fragments have lone pairs.".format(sum(use_lonepairs)),
                        italics=True)
    # Obtain the lists with the excluded pairs
    excluded_pairs_12 = config.getlist("fitting", "excluded_pairs_12", int)
    excluded_pairs_13 = config.getlist("fitting", "excluded_pairs_13", int)
    excluded_pairs_14 = config.getlist("fitting", "excluded_pairs_14", int)
    
    #Obtain charges (in the order of input), pols and polfacs
    charges = config.getlist("fitting", "charges", float)
    polarizabilities = config.getlist("fitting", "polarizabilities", float)
    polarizability_factors = config.getlist("fitting", "polarizability_factors", float)
    
    ################################################################################
    ## DIMER   PROPERTIES ##########################################################
    ################################################################################
    
    # Obtain A, C6 and d6 from user in the same order as the given pairs AA, AB ...
    # By default, if not specified in the config.ini file, these lists
    # will be set to a list of the same length, but with all coefficients 0.0
    C6 = config.getlist("fitting", "c6", float)
    
    # Initialize A and b to 0, for now
    try:
        A_buck = config.getlist("fitting", "A", float)
        b_buck = config.getlist("fitting", "d6", float)
        d6 = b_buck
    except:
        A_buck = [0.0] * len(C6)
        b_buck = [0.0] * len(C6)
        d6 = [0.0] * len(C6)
        if number_of_monomers < 3:
            system.format_print("You are using less than 3 monomers, but have no A or d6 constants in your input file.",
                                italics=True, color=system.Color.PURPLE)
            system.format_print("Are you sure this is correct?",
                                italics=True, color=system.Color.PURPLE)
    
    ################################################################################
    ## POLYNOMIAL PROPERTIES #######################################################
    ################################################################################
    
    #Ask the user for the max value of k and d
    # These are the non-linear parameters of the polynomials
    k_min = config.get("fitting", "k_min")
    k_max = config.get("fitting", "k_max")
    k_min_intra = k_min
    k_max_intra = k_max
    
    d_min = config.get("fitting", "d_min")
    d_max = config.get("fitting", "d_max")
    d_min_intra = d_min
    d_max_intra = d_max
    
    # Obtain inner and outer cutoff from config.ini
    ri = config.get("fitting", "r_in") # polynomials start to decrease (Angstrom)
    ro = config.get("fitting", "r_out") # polynomials start to decrease (Angstrom)
    
    # Define kind of variables for intra, inter and lone pairs
    # Options are:
    # exp e^-kd
    # exp0 e^-k(d-d0)
    # coul0 [e^-k(d-d0)]/r
    # coul [e^-kd]/r
    # Recomendation is to use exp for intra and inter and coul for lone pairs
    var_intra = config.get("fitting", "var_intra")
    var_virtual_sites = config.get("fitting", "var_virtual_sites")
    var_inter = config.get("fitting", "var_inter")
    
    system.format_print("Using {} variables for intramolecular interactions.".format(var_intra),
                        italics=True)
    system.format_print("Using {} variables for intermolecular interactions.".format(var_inter),
                        italics=True)
    system.format_print("Using {} variables for lone pair interactions.".format(var_virtual_sites),
                        italics=True)
    
    npoly = config.getint("fitting", "npoly")
    
    system.format_print("{} terms in the polynomial.".format(npoly),
                         italics=True)
    
    # Define Energy Range for the fitting
    E_range = 20.0
    
    # Define alpha for the fitting
    alpha = 0.0005
    
    ################################################################################
    ## Prepare pair information ####################################################
    ################################################################################
    
    # create dictionary mapping from atom index to atom name
    atom_list = [sym + str(ind) for sym, ind, frag in symmetry_parser.get_atoms()]
    
    system.format_print("Atoms in the molecule: {}.".format(" ".join(atom_list)),
                        italics=True)
    
    # read the poly.in file to get the variables
    variables, intra_poly_pairs, inter_poly_pairs = utils_nb_fitting.read_poly_in(name, virtual_sites_poly, var_intra, var_inter, var_virtual_sites)
    nvars = len(variables)
    
    system.format_print("{} variables in the polynomial.".format(nvars))
    
    ################################################################################
    ## Prepare non-linear terms information ########################################
    ################################################################################
    
    # Get the different non-linear parameters
    nlparam_intra, nlparam_inter, nl_param_all = utils_nb_fitting.get_non_linear_parameters(variables)
    # Save number in num_nonlinear
    num_nonlinear = len(nl_param_all)
    
    system.format_print("{} non-linear parameters in the polynomial.".format(num_nonlinear))
    
    
    ################################################################################
    ## Write monomer classes header and cpp ########################################
    ################################################################################
    
    # Write the header files
    for i in range(number_of_monomers):
        file_writer_nb_fitting.write_monomer_class_header(i + 1)
    
    # Write the cpp files
    for i in range(number_of_monomers):
        file_writer_nb_fitting.write_monomer_class_cpp(i + 1, number_of_sites[i], number_of_atoms[i], excluded_pairs_12[i],
                                                        excluded_pairs_13[i], excluded_pairs_14[i], charges[i],
                                                        polarizabilities[i], polarizability_factors[i])
    
    # Write mb-pol water monomer if applicable
    for i in range(number_of_monomers):
        if use_mbpol[i] == 1:
            file_writer_nb_fitting.write_mbpol_monomer(i + 1)
    
    
    ################################################################################
    ## Write polynomial holders header and cpp #####################################
    ################################################################################
    
    system_name = symmetry_parser.get_symmetry()
    
    file_writer_nb_fitting.write_fit_polynomial_holder_header(system_name, number_of_monomers, nl_param_all, ri, ro)
    
    file_writer_nb_fitting.write_fit_polynomial_holder_cpp(system_name, symmetry_parser, number_of_monomers, number_of_atoms, virtual_sites_poly, use_lonepairs, nl_param_all, variables, nvars, ri, ro, k_min_intra, k_max_intra, k_min, k_max, d_min_intra, d_max_intra, d_min, d_max)
    
    
    ################################################################################
    ## Write dispersion header and cpp #############################################
    ################################################################################
    if (number_of_monomers<3):
        file_writer_nb_fitting.write_dispersion_header(symmetry_parser, virtual_sites_poly, C6, d6)
        file_writer_nb_fitting.write_dispersion_cpp(symmetry_parser, virtual_sites_poly, excluded_pairs_12[0], excluded_pairs_13[0],excluded_pairs_14[0])
    
    
    
    ################################################################################
    ## Write buckingham header and cpp #############################################
    ################################################################################
        file_writer_nb_fitting.write_buckingham_header(symmetry_parser, virtual_sites_poly, A_buck, b_buck)
        file_writer_nb_fitting.write_buckingham_cpp(symmetry_parser, virtual_sites_poly, excluded_pairs_12[0], excluded_pairs_13[0],excluded_pairs_14[0])
    
    
    ################################################################################
    ## Polynomial header and cpp ###################################################
    ################################################################################
    
    file_writer_nb_fitting.write_poly_fit_header(number_of_monomers, system_name, degree, nvars, npoly)
    
    file_writer_nb_fitting.write_poly_fit_cpp(number_of_monomers, system_name, nvars, npoly, directcpp)
    
    ################################################################################
    ## Fitting code ################################################################
    ################################################################################
    
    file_writer_nb_fitting.write_fitting_code(number_of_monomers, number_of_atoms, number_of_sites, system_name, nl_param_all, k_min, k_max, d_min, d_max, k_min_intra, k_max_intra, d_min_intra, d_max_intra, E_range, alpha)
    
    ################################################################################
    ## Evaluation code #############################################################
    ################################################################################
    
    file_writer_nb_fitting.write_eval_code(number_of_monomers, number_of_atoms, number_of_sites, system_name)
    
    ################################################################################
    ## TTM-nrg fitting code ########################################################
    ################################################################################
    if number_of_monomers == 2:
        file_writer_nb_fitting.write_fitting_ttm_code(symmetry_parser, virtual_sites_poly, number_of_monomers, number_of_atoms, number_of_sites, system_name, k_min, k_max, E_range)
    
    ################################################################################
    ## TTM-nrg evaluation code #####################################################
    ################################################################################
    
        file_writer_nb_fitting.write_eval_ttm_code(symmetry_parser, virtual_sites_poly, number_of_monomers, number_of_atoms, number_of_sites, system_name)
    
    ################################################################################
    ## Makefile ####################################################################
    ################################################################################
    
    file_writer_nb_fitting.write_makefile(number_of_monomers, system_name)
    
    ################################################################################
    ## Write polynomial header and cpp for MBX #####################################
    ################################################################################
    
    file_writer_nb_fitting.write_poly_header_mbx(number_of_monomers, system_name, degree, nvars, npoly)
    
    if use_direct:
        file_writer_nb_fitting.write_direct_poly_cpp_grad_mbx(number_of_monomers, system_name, degree, nvars, npoly, poly_directory)

        file_writer_nb_fitting.write_direct_poly_cpp_nograd_mbx(number_of_monomers, system_name, degree, nvars, npoly, poly_directory)
    else:
        file_writer_nb_fitting.write_poly_cpp_grad_mbx(number_of_monomers, system_name, degree, nvars, npoly, poly_directory)
        
        file_writer_nb_fitting.write_poly_cpp_nograd_mbx(number_of_monomers, system_name, degree, nvars, npoly, poly_directory)
    
    ################################################################################
    ## Write polynomial holder and cpp for MBX #####################################
    ################################################################################
    
    file_writer_nb_fitting.write_mbx_polynomial_holder_header(number_of_monomers, system_name, degree, nvars, npoly, poly_directory, nl_param_all, ri, ro, virtual_sites_poly, "v1")
    
    file_writer_nb_fitting.write_mbx_polynomial_holder_cpp(system_name, symmetry_parser, number_of_monomers, number_of_atoms, virtual_sites_poly, use_lonepairs, nl_param_all, variables, nvars, degree, ri, ro, k_min_intra, k_max_intra, k_min, k_max, d_min_intra, d_max_intra, d_min, d_max, "v1")

