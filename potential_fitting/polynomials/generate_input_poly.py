# external package imports
import os

# absolute module imports
from potential_fitting.utils import SettingsReader, files, system
from potential_fitting.exceptions import InvalidValueError, InconsistentValueError

# local module imports
from .molecule_in_parser import MoleculeSymmetryParser

def generate_input_poly(settings_file, molecule_in, in_file_path, virtual_sites=["X", "Y", "Z"]):
    """
    Generates an input file for the polynomial generator.

    Args:
        settings_file       - Local path to the ".ini" file with all relevant settings.
        molecule_in         - String indicating the symmetry of the molecule, for example 'A1B2_A1B2'
        in_file_path       - Local path to the ".in" file to write the polynomial input file to.
        virtual_sites       - List of Symmetry labels that are virtual sites.
                Default: ["X", "Y", "Z"]
    """

    settings = SettingsReader(settings_file)

    system.format_print("Generating polynomial input file for symmetry {} into file {}.".format(molecule_in, in_file_path),
            bold=True, color=system.Color.YELLOW)

    # file will be automatically closed after block by with open as ... syntax
    with open(files.init_file(in_file_path, files.OverwriteMethod.NONE), "w") as poly_in_file:

        symmetry_parser = MoleculeSymmetryParser(molecule_in, virtual_sites=virtual_sites)

        # loop thru each fragment and add an add_molecule line at top of file
        for fragment_symmetry in symmetry_parser.get_fragment_symmetries():
            if "_" in fragment_symmetry:
                poly_in_file.write("add_molecule['({})']\n".format(fragment_symmetry))
            else:
                poly_in_file.write("add_molecule['{}']\n".format(fragment_symmetry))


        # newline between fragments and variables
        poly_in_file.write("\n")

        variable_count = 0

        for variable in symmetry_parser.get_variables():
            if not "intra" in variable[6] or not any([(vsite in variable[0] or vsite in variable[3]) for vsite in virtual_sites]):
                poly_in_file.write("add_variable['{}', '{}', '{}', '{}', '{}', '{}', '{}']\n".format(*variable))
                variable_count += 1

        # newline between variables and filters
        poly_in_file.write("\n")

        # add filter based on the filtering setting in settings.ini
        polynomial_filtering = settings.get("poly_generation", "accepted_terms",
                "all" if symmetry_parser.get_num_fragments() == 1 else "partly-inter")

        # make sure the user hasn't chosen intermolecular terms only with a monomer
        if symmetry_parser.get_num_fragments() == 1 and not polynomial_filtering == "all":
            raise InconsistentValueError("number of fragments", "[poly_generation].accepted_terms",
                    symmetry_parser.get_num_fragments(), polynomial_filtering,
                    "when there is only 1 fragment, there are no intermolecular interactions, "
                    + "so you must use 'all' terms")

        if polynomial_filtering == "purely-inter":
            # this filter filters out all terms that have any intra-molecular components
            poly_in_file.write("add_filter['ind-degree', 'x-intra-*+*-*', '1+']")
            system.format_print("Adding filter to filter out terms that use ANY intramolecular variables.",
                    italics=True)
        elif polynomial_filtering == "partly-inter":
            # this filter filters out all terms that have no inter-molecular components
            poly_in_file.write("add_filter['sum-degree', 'x-inter-*+*-*', '0']")
            system.format_print("Adding filter to filter out terms that ONLY use intramolecular variables.",
                    italics=True)

        elif polynomial_filtering != "all":
            raise InvalidValueError("[poly_generation][accepted_terms]",polynomial_filtering,
                    "one of 'all', 'partly-inter', or 'purely-inter'")

    system.format_print("Successfully generated polynomial input file! {} total variables.".format(variable_count),
            bold=True, color=system.Color.GREEN)
