import sys, os, json
from potential_fitting.utils import SettingsReader
from potential_fitting.utils import constants, system
from potential_fitting.polynomials import MoleculeSymmetryParser
from . import utils_nb_fitting 
from . import file_writer_nb_fitting 


def generate_ttmnrg_fitting_code(settings_path, config_path):

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
    
    A_buck = [0.0] * len(C6)
    b_buck = [0.0] * len(C6)
    d6 = [0.0] * len(C6)
    
    ################################################################################
    ## POLYNOMIAL PROPERTIES #######################################################
    ################################################################################
    
    #Ask the user for the max value of k and d
    # These are the non-linear parameters of the polynomials
    d_min = config.get("fitting", "d_min")
    d_max = config.get("fitting", "d_max")
    d_min_intra = d_min
    d_max_intra = d_max
    
    # Define Energy Range for the fitting
    E_range = config.get("fitting","energy_range",float)

    # Define alpha for the fitting
    alpha = config.get("fitting","alpha",float)
    
    ################################################################################
    ## Prepare pair information ####################################################
    ################################################################################
    
    # create dictionary mapping from atom index to atom name
    atom_list = [sym + str(ind) for sym, ind, frag in symmetry_parser.get_atoms()]
    
    system.format_print("Atoms in the molecule: {}.".format(" ".join(atom_list)),
                        italics=True)
    
    
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
    
    ################################################################################
    ## Write dispersion header and cpp #############################################
    ################################################################################
    file_writer_nb_fitting.write_dispersion_header(symmetry_parser, virtual_sites_poly, C6, d6)
    file_writer_nb_fitting.write_dispersion_cpp(symmetry_parser, virtual_sites_poly, excluded_pairs_12[0], excluded_pairs_13[0],excluded_pairs_14[0])
    
    
    
    ################################################################################
    ## Write buckingham header and cpp #############################################
    ################################################################################
    file_writer_nb_fitting.write_buckingham_header(symmetry_parser, virtual_sites_poly, A_buck, b_buck)
    file_writer_nb_fitting.write_buckingham_cpp(symmetry_parser, virtual_sites_poly, excluded_pairs_12[0], excluded_pairs_13[0],excluded_pairs_14[0])
   
    
    ################################################################################
    ## TTM-nrg fitting code ########################################################
    ################################################################################

    file_writer_nb_fitting.write_fitting_ttm_code(symmetry_parser, virtual_sites_poly, number_of_monomers, number_of_atoms, number_of_sites, system_name, k_min, k_max, E_range)
    
    ################################################################################
    ## TTM-nrg evaluation code #####################################################
    ################################################################################
    
    file_writer_nb_fitting.write_eval_ttm_code(symmetry_parser, virtual_sites_poly, number_of_monomers, number_of_atoms, number_of_sites, system_name)
    
    ################################################################################
    ## Makefile ####################################################################
    ################################################################################
    
    file_writer_nb_fitting.write_makefile_ttmnrg(number_of_monomers, system_name)
    
