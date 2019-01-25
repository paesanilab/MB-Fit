# external package imports
import os

# absolute module imports
from potential_fitting.utils import SettingsReader, files
from potential_fitting.exceptions import InvalidValueError, InconsistentValueError

# local module imports
from .molecule_in_parser import MoleculeInParser

def generate_input_poly(settings_file, molecule_in, in_file_path):
    """
    Generates an input file for the polynomial generator.

    Args:
        settings_file       - Local path to the ".ini" file with all relevant settings.
        molecule_in         - String indicating the symmetry of the molecule, for example 'A1B2_A1B2'
        in_file_path       - Local path to the ".in" file to write the polynomial input file to.
    """

    settings = SettingsReader(settings_file)

    # file will be automatically closed after block by with open as ... syntax
    with open(files.init_file(in_file_path, files.OverwriteMethod.NONE), "w") as poly_in_file:

        molecule_in_parser = MoleculeInParser(molecule_in)

        # loop thru each fragment and add an add_molecule line at top of file
        for fragment_parser in molecule_in_parser.get_fragments():
            poly_in_file.write("add_molecule['{}']\n".format(fragment_parser.get_fragment_in()))

        # newline between fragments and variables
        poly_in_file.write("\n")

        # these letters mark virtual sites instead of atoms
        virtual_sites = ['X', 'Y', 'Z']

        for variable in molecule_in_parser.get_variables():
            poly_in_file.write("add_variable['{}', '{}', '{}', '{}', '{}']\n".format(*variable))

        # newline between variables and filters
        poly_in_file.write("\n")

        # add filter based on the filtering setting in settings.ini
        polynomial_filtering = settings.get("poly_generation", "accepted_terms",
                "all" if len(list(molecule_in_parser.get_fragments())) == 1 else "purely-inter")

        # make sure the user hasn't chosed intermolecular terms only with a monomer
        if len(list(molecule_in_parser.get_fragments())) == 1 and not polynomial_filtering == "all":
            raise InconsistentValueError("number of fragments", "[poly_generation].accepted_terms",
                    len(list(molecule_in_parser.get_fragments())), polynomial_filtering,
                    "when there is only 1 fragment, there are no intermolecular interactions, "
                    + "so you must use 'all' terms")

        if polynomial_filtering == "purely-inter":
            # this filter filters out all terms that have any intra-molecular components
            poly_in_file.write("add_filter['degree', 'x-intra-*+*', '1+', '*']")
        elif polynomial_filtering == "partly-inter":
            # this filter filters out all terms that have no inter-molecular components
            poly_in_file.write("add_filter['not', 'degree', 'x-inter-*+*', '1+', '*']")

        elif polynomial_filtering != "all":
            raise InvalidValueError("[poly_generation][accepted_terms]",polynomial_filtering,
                    "one of 'all', 'partly-inter', or 'purely-inter'")
