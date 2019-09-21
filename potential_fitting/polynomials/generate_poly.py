# external package imports
import itertools

# absolute module imports
from potential_fitting.utils import SettingsReader, files, system, ProgressBar
from potential_fitting.exceptions import ParsingError, InvalidInputError, InconsistentValueError

# relative module imports
from . import filters
from .molecule_in_parser import MoleculeInParser, FragmentParser
from .monomial import Monomial
from .variable import Variable


class PolynomialGenerator(object):
    """
    Class for generating of polynomial monomials and cpp and maple files.
    """

    def __init__(self, settings_path, num_gradient_terms_per_line=1):
        """
        Constructs a new PolynomialGenerator.

        Args:
            settings_path   - Local path to '.ini' settings file containing settings for this PolynomialGenerator.
            num_gradient_terms_per_line - The number of terms to put in each line in the cpp file
                    for computing the gradients. Larger numbers may result in slower compilation times.

        Returns:
            A new PolynomialGenerator.
        """

        self.settings = SettingsReader(settings_path)
        self.num_gradient_terms_per_line = num_gradient_terms_per_line

    def generate_polynomial(self, input_path, order, output_dir_path):

        """
        Generates .h, .cpp, .maple, and .log polynomial files for a new polynomial based on the input file and order.

        Args:
            input_path      - Local path to '.in' polynomial input file.
            order           - Order of the polynomial to generate.
            output_dir_path - Local path to the directory in which to write all files.

        Returns:
            None.
        """

        system.format_print("Generating polynomial files from input file {} into directory {}.".format(input_path, output_dir_path),
                            bold=True, color=system.Color.YELLOW)

        log_path = "{}/poly.log".format(output_dir_path)

        monomials, variables, variable_permutations = self.get_monomials_and_variables(input_path, order, log_path)

        self.write_cpp_files(monomials, variables, variable_permutations, output_dir_path)

        system.format_print("Successfully generated polynomial files into directory {}.".format(output_dir_path),
                            bold=True, color=system.Color.GREEN)

    def get_monomials_and_variables(self, input_path, order, log_path):
        """
        Gets Monomials, variables, and variable permutations for a polynomial defined by the input file and given order.

        Args:
            input_path      - Local path to '.in' polynomial input file.
            order           - Order of monomials to generate.
            log_path        - Local path to '.log' file to write the log to.

        Returns:
            (monomials, variables, permutations)
            monomials       - List of all Monomials in the polynomial.
            variables       - List of all variables in the polynomial.
            permutations    - List of all the permutations of the variables in this polynomial that can be used to
                    permute the monomials to get the permutations of each Monomial.
        """

        system.format_print("Parsing polynomial input file {}...".format(input_path),
                            italics=True)

        # parse fragments, variables, and filters from the input file
        fragments, variables, monomial_filters = self.read_input_file(input_path)

        # create a molecule parser of all the fragments.
        molecule_in_parser = MoleculeInParser("_".join(fragments))

        # Parse the fragment names to name the atoms in the system

        # list of names of atoms in the molecule, each name consists of a type (ex: 'A', 'B', etc) concatenated to
        # a letter that indicates the fragment (ex: 'a', 'b', etc).
        atom_names = []

        # holds the index of the first atom in each fragment.
        index_each_fragment = []

        # index of the next atom to parse.
        index = 0

        for fragment_parser in molecule_in_parser.get_fragments():
            index_each_fragment.append(index)
            for atom_type_parser in fragment_parser.get_atom_and_virtual_site_types():
                for atom in atom_type_parser.get_atoms():
                    atom_names.append("{}{}".format(atom, fragment_parser.get_fragment_id()))
                    index += 1

        print("Atom names", atom_names)

        files.init_file(log_path)

        # open the log file
        with open(log_path, "w") as log_file:

            # write number of fragments to log
            log_file.write("<> fragments({}) <>\n".format(len(fragments)))
            log_file.write("\n")

            # write each fragment to log
            for fragment in fragments:
                log_file.write(fragment + "\n")
            log_file.write("\n")

            # write number of atoms to log
            log_file.write("<> atoms ({}) <>\n".format(len(atom_names)))
            log_file.write("\n")

            # write each atom to log
            log_file.write(":".join(atom_names) + "\n")
            log_file.write("\n")

            system.format_print("Finding permutations...",
                                italics=True)

            atom_permutations = list(self.get_permutations("_".join(fragments)))

            # write permutation count to log file
            log_file.write("<> permutations ({}) <>\n".format(len(atom_permutations)))
            log_file.write("\n")

            # write each permutation to log file
            for permutation in atom_permutations:
                log_file.write(":".join(str(x) for x in permutation) + "\n")
            log_file.write("\n")

            # write variable count to log file
            log_file.write("<> variables ({}) <>\n".format(len(variables)))
            log_file.write("\n")

            # generate the variable permutations
            variable_permutations = list(self.get_variable_permutations(variables, atom_permutations, atom_names))

            # write each variable to log file
            for index, variable in zip(range(len(variables)), variables):
                log_file.write("{:>2} : {:>3}({:>1}) <===> {:>3}({:>1}) : {}\n".format(index, variable.atom1_name,
                                                                                            variable.atom1_fragment,
                                                                                            variable.atom2_name,
                                                                                            variable.atom2_fragment,
                                                                                            variable.category))
            log_file.write("\n")

            # generate the monomials

            # this list will contain all accepted monomials for the polynomial.
            total_monomials = []

            system.format_print("Generating terms up to degree {}...".format(order),
                                italics=True)

            unfiltered_monomials = self.get_monomials_dynamic(len(variables), order, variable_permutations)

            # loop thru every degree in this polynomial
            for degree in range(1, order + 1):

                # header for this degree
                log_file.write("<> {} degree <>\n".format(degree))
                log_file.write("\n")

                # get all the monomials of the current degree
                monomials = unfiltered_monomials[degree - 1]

                # log number of possible monomials
                log_file.write("{} possible {} degree monomials\n".format(len(monomials), degree))

                system.format_print("{} possible degree {} terms, now filtering them...".format(len(monomials), degree),
                                    italics=True)

                # filter out monomials from the list based on filters in the poly.in file
                accepted_monomials = list(self.filter_monomials(monomials, variables, monomial_filters))

                # log number of accpeted terms
                log_file.write("{} <<== accepted {} degree terms\n".format(len(accepted_monomials), degree))

                log_file.write("\n")

                # add the monomials of the current degree to the list of all monomials
                total_monomials.extend(accepted_monomials)

                system.format_print("There were {} accepted degree {} terms.".format(len(accepted_monomials), degree),
                                    italics=True)

            system.format_print("There were {} accepted terms over all".format(len(total_monomials)),
                                italics=True)

            # log the total number of terms
            log_file.write(" Total number of terms: {}\n".format(len(total_monomials)))

            return total_monomials, variables, variable_permutations

    def write_cpp_files(self, monomials, variables, variable_permutations, output_dir_path):
        """
        Write the .cpp, .h, and .maple files for a polynomial.

        Args:
            monomials       - List of all Monomials in the polynomial.
            variables       - List of all the Variables in the polynomial.
            variable_permutations - List of all permutations of the Variables in this polynomial that can be used to
                    permute the monomials to get the permutations of each Monomial.
            output_dir_path - Local path to directory in which to write output files.

        Returns:
            None.
        """

        system.format_print("Writing .h and .maple polynomial files in directory {}...".format(output_dir_path),
                            italics=True)

        files.init_directory(output_dir_path)

        # make the .cpp vars file
        with open(output_dir_path + "/vars.cpp", "w") as vars_file:
            self.write_variable_file(vars_file, variables)

        # write the header file
        with open(output_dir_path + "/poly-model.h", "w") as header_file:
            self.write_header_file(header_file, len(monomials), len(variables))

        # open the three files we will now write to
        with open(output_dir_path + "/poly-direct.cpp", "w") as cpp_file, \
             open(output_dir_path + "/poly-grd-direct.cpp", "w") as cpp_grd_file, \
             open(output_dir_path + "/poly-grd.maple", "w") as grd_file, \
             open(output_dir_path + "/poly-nogrd.maple", "w") as nogrd_file:

            # write the opening for the cpp file
            self.write_cpp_opening(cpp_file, len(monomials), len(variables))
            self.write_cpp_grd_opening(cpp_grd_file, len(monomials), len(variables))

            system.format_print('Writing gradients...',
                                italics=True)
            progress_bar = ProgressBar(start=0, end=len(variables) - 1)
            progress_bar.update(0)

            for index in range(len(variables)):
                self.write_cpp_gradient(cpp_grd_file, index, monomials, variable_permutations)
                progress_bar.update(index)
            progress_bar.finish()

            # keeps track of what index in a list of all monomials the current monomial would occupy
            monomial_index = 0

            system.format_print('Writing terms...',
                                italics=True)
            progress_bar = ProgressBar(start=0, end=len(monomials))
            progress_bar.update(0)

            # loop thru every degree in this polynomial
            for monomial in monomials:
                self.write_cpp_monomial(cpp_file, monomial_index, monomial, variable_permutations)
                self.write_cpp_monomial(cpp_grd_file, monomial_index, monomial, variable_permutations)
                self.write_grd_monomial(grd_file, monomial_index, monomial, variable_permutations)
                self.write_nogrd_monomial(nogrd_file, monomial_index, monomial, variable_permutations)
                monomial_index += 1
                progress_bar.progress(1)

            progress_bar.finish()

            cpp_file.write("\n")
            grd_file.write("\n")
            nogrd_file.write("\n")

            self.write_cpp_closing(cpp_file, len(monomials))
            self.write_cpp_grd_closing(cpp_grd_file, len(monomials))
            self.write_grd_closing(grd_file, len(monomials), len(variables))
            self.write_nogrd_closing(nogrd_file, len(monomials), len(variables))

        system.format_print("Successfully generated polynomial files!",
                            italics=True)

    def read_input_file(self, input_path):
        """
        Reads the polynomial input file and returns a list of fragments, Variables, and Filters.

        Args:
            input_path          - Local path to the ".in" polynomial input file.

        Returns:
            A 3-tuple (fragments, variables, filters) as read from the input file. fragments is a list of "A1B2" type
            formats, variables is a list of Variable objects, filters is a list of Fitler objects.
        """

        fragments = []
        variables = []
        monomial_filters = []

        with open(input_path, "r") as input_file:

            # loop thru all lines in the input file
            for line in input_file:

                # if line starts with a pound sign, it is a comment line, so ignore it.
                if line.startswith("#"):
                    continue

                # check if this line is an add_molecule statement
                if line.startswith("add_molecule"):

                    # parse line into a fragment
                    fragments.append(line[line.index("['") + 2:line.index("']")])
                    continue

                # check if this line is an add_variable statement
                if line.startswith("add_variable"):

                    # remove all the extra information before and after the brackets
                    line = line[line.index("[") + 1:line.index("]")]

                    # remove all single quotes and spaces from  line
                    line = line.replace("'", "").replace(" ", "")

                    # split line into arguments and construct a new variable
                    variables.append(Variable(*line.split(",")))
                    continue

                # check if this line is an add_filter statement
                if line.startswith("add_filter"):

                    # remove all the extra information before and after the brackets
                    line = line[line.index("[") + 1:line.index("]")]

                    # remove all single quotes and spaces from  line
                    line = line.replace("'", "").replace(" ", "")

                    # split line into arguments and construct a new filter
                    monomial_filters.append(filters.parse_filter(*line.split(",")))
                    continue

                # unless this line is blank, it is an invalid format
                if line != "\n":
                    raise ParsingError(input_path, "Unrecongnized line format: '{}'".format(line))

        return fragments, variables, monomial_filters

    def get_permutations(self, molecule_in):

        if "_" in molecule_in or "(" in molecule_in or ")" in molecule_in:
            fragment_ins = self.split_molecule_in(molecule_in)
            fragment_permutations_list = [list(self.get_permutations(fragment_in)) for fragment_in in fragment_ins]
            fragment_atom_index = 0
            for fragment_permutations in fragment_permutations_list:
                for permutation in fragment_permutations:
                    for i in range(len(permutation)):
                        permutation[i] += fragment_atom_index
                fragment_atom_index += len(fragment_permutations[0])

            return self.get_molecule_permutations(fragment_ins, fragment_permutations_list)
        else:
            parser = FragmentParser(molecule_in, 'a')
            return self.get_fragment_permutations(parser.get_atom_and_virtual_site_types())
        # TODO FINISH THIS

    def split_molecule_in(self, molecule_in):
        fragments = []
        cur_fragment = ""

        num_open_parens = 0

        for char in molecule_in:

            if char is '(':
                if num_open_parens > 0:
                    cur_fragment += char
                num_open_parens += 1

            elif char is ')':
                num_open_parens -= 1
                if num_open_parens > 0:
                    cur_fragment += char

            elif char is '_' and num_open_parens == 0:
                fragments.append(cur_fragment[:])
                cur_fragment = ""

            else:
                cur_fragment += char

        if num_open_parens != 0:
            raise InvalidInputError("{} is not a valid symmetry input. Make sure all parenthesis are closed!".format(molecule_in))

        fragments.append(cur_fragment)
        return fragments

    def get_fragment_permutations(self, atom_type_generator):
        """
        Takes in a fragment and constructs all permutations of that fragment by swapping positions of equivalent atoms.

        Example Input and output:
        A3 -> [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]
        A1B2 -> [[0, 1, 2], [0, 2, 1]
        A1B1C1 -> [[0, 1, 2]]
        A2B2 -> [[0, 1, 2, 3], [0, 1, 3, 2], [1, 0, 2, 3], [1, 0, 3, 2]]
        A1B1A1 -> [[0, 1, 2],[2, 1, 0]]

        Args:
            atom_types            - Generator that generates one AtomParser for each atom type in the fragment.
                    Like those created by FragmentParser.get_atom_and_virtual_site_types().

        Generates:
            Lists of indices representing permutations of the fragment, created by forming all permutations of
            equivalent atoms. If [0, 1] and [1, 0] are both included in the output, then atoms 0 and 1 are interchangeable.
        """

        try:
            # get the first atom type from the generator.
            atom_type_parser = next(atom_type_generator)
        except StopIteration:
            # If the generator is exhausted, then there are no atoms in this fragment, so there are no permutations.
            yield []
            return

        # how many atoms of the first type are in the fragment?
        count = atom_type_parser.get_count()

        # get all possible permutations of the atoms of the first atom type.
        first_permutations = itertools.permutations(range(count))

        # convert first_permutations to a list of lists
        first_permutations = [list(first_permutation) for first_permutation in first_permutations]

        # get all possible permutations of the atoms of all atom types except the first.
        other_permutations = self.get_fragment_permutations(atom_type_generator)

        # offset all other_permutations by the number of atoms of the first type, also convert other_permutations to list of lists.
        other_permutations = [[count + i for i in perm] for perm in other_permutations]

        # yield every combination of permutations of the atoms of the first type and permuations of the atoms of all other types.
        yield from (first_permutation + other_permutation
                    for first_permutation, other_permutation
                    in itertools.product(first_permutations, other_permutations))

    def get_molecule_permutations(self, fragments, fragment_permutations):
        """
        Takes a list of fragments ("A1B2", "A1B1C1", etc) and a list of the permutations of those fragments and yields
        permutations of the molecule as a whole.

        Works by generating one permutation for every possible choice of 1 permutation from each fragment, but if 2
        fragments are the same, their permutations will also be swapped with each other to create all permutations where
        equivalent atoms are swapped in every possible order.

        Args:
            fragments           - List of strings of format "A1B2" representing the fragments.
            fragment_permutations - List of permutations for each fragment as generated by get_fragment_permutations().

        Generates:
            Lists of indices representing all permutations of the atoms in this molecule where equivalent atoms
            are swapped in every possible order.
        """

        if len(fragments) != len(fragment_permutations):
            raise InconsistentValueError("number of fragments", "number of fragment permutations",
                                         len(fragments), len(fragment_permutations),
                                         "these values must be equal")

        # if there are no fragments, then there are no permutations for this molecule
        if len(fragments) == 0:
            yield []
            return

        # get list of permutations of the first fragment
        first_permutations = fragment_permutations[0]

        # get list of permutations of every other fragment in the molecule
        other_permutations = self.get_molecule_permutations(fragments[1:], fragment_permutations[1:])

        # yield every combination of permutations of the first fragment and permutations of the rest of the molecule.
        yield from (first_permutation + other_permutation
                    for first_permutation, other_permutation
                    in itertools.product(first_permutations, other_permutations))

        # loop thru every other fragment in the molecule
        for i in range(1, len(fragments)):

            # check if the first fragment is the same the other fragment
            # if so, we must yield additional permutations where these fragments are switched.
            if fragments[i] == fragments[0]:

                # get the list of all permutations of the other fragment that is the same as the first one.
                other_permutations = fragment_permutations[i]

                # get list of permutations of every fragment except the other fragment that is the same as the first one
                # after swapping the first fragment with the other one.
                switched_permutations = self.get_molecule_permutations(fragments[1:i] + [fragments[0]] + fragments[i + 1:],
                                                             fragment_permutations[1:i] + [fragment_permutations[0]] +
                                                             fragment_permutations[i + 1:])

                # yield every combination of permutations of the other fragment and the permutations produced by
                # switching the first fragment with the other one.
                yield from (other_permutation + switched_permutation
                            for other_permutation, switched_permutation
                            in itertools.product(other_permutations, switched_permutations))

    def get_variable_permutations(self, variables, molecule_permutations, atom_names):
        """
        Constructs the variable permutations from the molecule permutations.

        Args:
            variables           - A list of the variables in this polynomial.
            molecule_permutations   - A list of all the molecule permutations as generated by get_molecule_permutations().
            atom_names          - A list of the names of each atom in this molecule, unique for each atom. Must match
                    those in the variables.

        Generates:
            Lists of indices representing the variable permutations in this fragment. Where if permutation[0] = 1, then
            in that permutation, variable at index 1 is moved to index 0.
        """

        # there will be 1 variable permutation for every molecule permutation
        for atom_permutation in molecule_permutations:

            # initialize the new variable permutation
            variable_permutation = []

            # each variable permutation has length == len(variables) where each value is a value of a variable which should
            # be switched with this index to create this permutation
            for variable in variables:

                # get the variable's atoms
                atom1 = variable.atom1_name + variable.atom1_fragment
                atom2 = variable.atom2_name + variable.atom2_fragment

                # get the permutated variable's atoms

                new_atom1 = atom_names[atom_permutation[atom_names.index(atom1)]]
                new_atom2 = atom_names[atom_permutation[atom_names.index(atom2)]]

                # find the index of this variable in the list of variables
                new_index = -1
                for new_variable_index, new_variable in enumerate(variables):
                    # if both new atoms match the atoms for a variable, then that is the permutated variable
                    if (new_atom1 == new_variable.atom1_name + new_variable.atom1_fragment and new_atom2 ==
                            new_variable.atom2_name + new_variable.atom2_fragment):
                        new_index = new_variable_index
                        break
                    if (new_atom2 == new_variable.atom1_name + new_variable.atom1_fragment and new_atom1 ==
                            new_variable.atom2_name + new_variable.atom2_fragment):
                        new_index = new_variable_index
                        break

                # this catch is technically not needed, as this should *never* happen.
                # a variable should always be found, if one is not found new_index being -1 signifies a problem
                if new_index == -1:
                    print("Uh Oh, something went wrong. This should never print. :( Contant Ethan to help debug this very cryptic error message.")
                    print("Its not a you problem, its a Ethan problem.")
                    raise Exception

                # append the index of the permutated variable to variable permutation
                variable_permutation.append(new_index)

            yield variable_permutation

    def get_monomials_dynamic(self, number_of_variables, degree, variable_permutations):

        # define monomial grid as 2d array. where x[n][d] means all monomials using only the first n
        # variables of total degree d.
        # each element x[n][d] will be a set of such monomials.
        monomial_grid = [[None for i in range(0, degree + 1)] for i in range(0, number_of_variables + 1)]

        # first we fill in our base cases

        # any monomial with total degree 0 has every variable raised to the 0th power
        d = 0
        for n in range(0, number_of_variables + 1):
            monomial_grid[n][d] = set([Monomial([0 for k in range(n)])])

        # There are no monomials with no variables and positve total degree, so all sets where n = 0 are empty.
        # (except the one where d = 0 as well, where we have one monomial with no variables)
        n = 0
        for d in range(1, degree + 1):
            monomial_grid[n][d] = set([])

        progress_bar = ProgressBar(start=0, end=(number_of_variables) * (degree))
        progress_bar.update(0)

        # Now, for our recursive case!

        # Fill in monomial_grid starting from n = 0, and ending at n = number_of_variables.
        # (Remember, the row n=0 was filled in as a base case.)
        for n in range(1, number_of_variables + 1):

            # get a list of all variable permutations using only the first n variables.
            new_var_permutations = [p[:n] for p in variable_permutations if all([i < n for i in p[:n]])]

            # Fill in each row of monomial_grid starting from d = 1 and ending at d = degree.
            # (Remember, the row d=0 was filled in as a base case.)
            for d in range(1, degree + 1):

                # initialize a set of permutationally independent monomials.
                monomials = set()

                # each monomial could have degree [0, d] in the last variable of the first n.
                # we loop over these degrees and for each case add monomials to the list of monomials.
                for i in range(0, d + 1):

                    # add all monomials with degree i in the last variable of the first n to the
                    # list of monomials.
                    # this step is done by looking at all monomials with d-i total degree and one less variable and
                    # appending an extra term of degree i to the end.
                    for monomial in set([Monomial(monomial.get_degrees() + [i])
                                        for monomial
                                        in monomial_grid[n-1][d-i]]):
                        
                        monomials.add(monomial.get_standard_permutations(new_var_permutations))

                # update the element of monomial_grid to equal all permutationally independent monomials of total
                # degree d using just the first n variables.
                monomial_grid[n][d] = monomials

                progress_bar.update((n-1)*degree + d)

        progress_bar.finish()

        # return a list of length degree L, where each element L[d] is a list of all monomials using all variables
        # with total degree d + 1. List will be of length degree.
        return [monomial_grid[number_of_variables][d] for d in range(1, degree + 1)]

    def get_monomials(self, number_of_variables, degree, variable_permutations):
        """
        Given a number of variables and a degree, generates all possible monomials of that degree.

        Yields every possible list of length number_of_vars with entries adding to degree.
        All entries are non-negative integers.

        Args:
            number_of_variables - The number of variables to generate monomials for.
            degree              - The degree of each monomial.

        Generates:
            All possible lists of length number of vars with entries adding to degree where each entry is a
            non-negative integer.
        """

        # if number of vars is 1, then the only possible monomial is a monomial with 1 variable with the given degree
        if number_of_variables == 1:
            yield Monomial([degree])
            return

        # if degree is 0, then the only possible monomial is a monomial with all 0 degrees with the given number of
        # variables
        if degree == 0:
            yield Monomial([0 for i in range(number_of_variables)])
            return

        yielded_monomials = set()

        # loop over all possible first non-zero terms
        for first_non_zero_term_index in range(number_of_variables):

            new_variable_permutations = [[i - first_non_zero_term_index - 1 for i in p[first_non_zero_term_index + 1:]]
                                        for p
                                        in variable_permutations
                                        if all([i - first_non_zero_term_index - 1 >= 0 for i in p[first_non_zero_term_index + 1:]])]

            # loop over all possible values of the first non-zero term's degree
            for first_non_zero_degree in range(1, degree + 1):

                # generate the list of degrees for all terms before the first non-zero term.
                zero_terms = [0 for i in range(first_non_zero_term_index)]

                # generate the list of degrees for the first non-zero term.
                first_non_zero_term = [first_non_zero_degree]

                # get all the possible lists of degrees for all terms after the first non-zero term.
                other_terms = self.get_monomials(number_of_variables - first_non_zero_term_index - 1,
                                                 degree - first_non_zero_degree,
                                                 new_variable_permutations)

                # yield every monomial created by concatenating all the terms before the first non-zero term, the first
                # non-zero term, and every possible list of terms after the first non-zero term.

                mons = (Monomial(zero_terms + first_non_zero_term + other_term.degrees)
                            for other_term
                            in other_terms)

                for mon in mons:

                    standard_mon = mon.get_standard_permutations(variable_permutations)

                    if standard_mon in yielded_monomials:
                        continue

                    yielded_monomials.add(standard_mon)
                    yield standard_mon

    def filter_monomials(self, monomials, variables, filters):
        """
        Filters a list of monomials with the given filters.

        Args:
            monomials           - List of Monomials as generated by get_monomials().
            variables           - List of all variables in the monomials.
            monomial_filters    - List of Filters to use to filter the monomials.

        Returns:
            List of Monomials, without any Monomial that fails at least 1 filter.
        """

        # we must try each monomial
        for monomial in monomials:

            keep = True

            # if the monomial fails 1 ore more filters, then we do not yield it.
            for filter in filters:
                if not filter.keep(monomial, variables):
                    keep = False
                    break

            # if the monomial didn't fail any filters, yield it.
            if keep:
                yield monomial

    def write_variable_file(self, variable_file, variables):
        """
        Writes the 'variable.cpp' file.

        Args:
            variable_file   - Local path to file to write in.
            variables       - All variables to write in the file.

        Returns:
            None.
        """
        # variable header comment
        variable_file.write("\n// <-> variables <->\n\n")

        # loop thru each variable
        for index, variable in enumerate(variables):

            # write the variable to vars.cpp
            variable_file.write("    x[{}] = @VAR@-|{}|({}{}{}{};\n".format(index, variable.category, variable.atom1_name,
                    variable.atom1_fragment, variable.atom2_name, variable.atom2_fragment))

    def write_header_file(self, header_file, total_terms, number_of_variables):
        """
        Writes the '.h' polynomial file.

        Args:
            header_file     - Local path to file to write in.
            total_terms     - Number of terms in the polynomial.
            number_of_variables - Number of variables in the polynomial.

        Returns:
            None.
        """
        header_file.write("""#ifndef POLY_MODEL_H
#define POLY_MODEL_H

namespace mb_system {{

struct poly_model {{
    static const unsigned n_vars = {1};
    static const unsigned size = {0};

    static double eval(const double a[{0}],
                       const double x[{1}]);

    static double eval(const double a[{0}],
                       const double x[{1}],
                             double g[{1}]);

    static double eval_direct(const double a[{0}],
                              const double x[{1}]);
    static double eval_direct(const double a[{0}],
                              const double x[{1}], 
                              double g[{1}]);

public:
    unsigned report_nvars(){{ return n_vars; }};
    unsigned report_size(){{ return size; }};
}};

}} // namespace mb_system

#endif // POLY_MODEL_H""".format(total_terms, number_of_variables))

    def write_cpp_opening(self, cpp_file, total_terms, number_of_variables):
        """
        Writes the opening of the 'poly-direct.cpp' polynomial file.

        Args:
            cpp_file        - Local path to file to write in.
            total_terms     - Number of terms in the polynomial.
            number_of_variables - Number of variables in the polynomial.

        Returns:
            None.
        """
        cpp_file.write("""#include "poly-model.h"
    
namespace mb_system {{

double poly_model::eval_direct(const double a[{0}], const double x[{1}])
{{
    double p[{0}];
""".format(total_terms, number_of_variables))

    def write_cpp_grd_opening(self, cpp_file, total_terms, number_of_variables):
        """
        Writes the opening of the 'poly-grd-direct.cpp' polynomial file.

        Args:
            cpp_file        - Local path to file to write in.
            total_terms     - Number of terms in the polynomial.
            number_of_variables - Number of variables in the polynomial.

        Returns:
            None.
        """
        cpp_file.write("""#include "poly-model.h"
    
namespace mb_system {{

double poly_model::eval_direct(const double a[{0}], const double x[{1}], double g[{1}])
{{
    double p[{0}];
""".format(total_terms, number_of_variables))

    def write_cpp_monomial(self, cpp_file, index, monomial, variable_permutations):
        """
        Writes a single monomial to the 'poly-direct.cpp' polynomial file.

        Args:
            cpp_file        - Local path to file to write in.
            index           - The index of this monomial.
            monomial        - The Monomial to write.
            variable_permutations - All permutations of the variables in this polynomial that can be used to permute
                    the Monomial.

        Returns:
            None.
        """
        cpp_file.write("    p[{}] = ".format(index))

        first_term = True
        for permutation in set(
                [permutation for permutation in monomial.permute(variable_permutations)]):

            if not first_term:
                cpp_file.write(" + ")

            first_term = False

            first_factor = True

            for variable_index, degree in enumerate(permutation):

                if degree == 0:
                    continue

                for term in range(degree):

                    if not first_factor:
                        cpp_file.write("*")

                    first_factor = False

                    cpp_file.write("x[{}]".format(variable_index))

        cpp_file.write(";\n")

    def write_cpp_gradient(self, cpp_file, index, monomials, variable_permutations):

        num_monomials = len(monomials)

        monomial_strings = []

        num_gradient_terms_on_line = 0
        t_index = 0

        for monomial_index, monomial in enumerate(monomials):

            permutation_strings = []

            for permutation in set(
                    [permutation for permutation in monomial.permute(variable_permutations)]):

                const, derivative = permutation.get_derivative(index)
                if const == 0:
                    continue

                variable_strings = [str(const)]

                for variable_index, degree in enumerate(derivative.degrees):
                    for i in range(degree):
                        variable_strings.append("x[{}]".format(variable_index))

                permutation_strings.append("*".join(variable_strings))

            if len(permutation_strings) != 0:
                monomial_strings.append("a[{}]*({})".format(monomial_index, " + ".join(permutation_strings)))

                num_gradient_terms_on_line += 1

                if num_gradient_terms_on_line == self.num_gradient_terms_per_line:

                    cpp_file.write("    const double t{}_{} = {};\n".format(index, t_index, " + ".join(monomial_strings)))

                    monomial_strings = []
                    num_gradient_terms_on_line = 0

                    t_index += 1

        if num_gradient_terms_on_line > 0:

            cpp_file.write("    const double t{}_{} = {};\n".format(index, t_index, " + ".join(monomial_strings)))

            t_index += 1

        end_t_index = t_index

        if end_t_index == 0:
            gradient_string = "    g[{}] = 0;\n".format(index)
        else:
            gradient_string = "    g[{}] = {};\n".format(index, " + ".join(["t{}_{}".format(index, t_index) for t_index in range(0, end_t_index)]))

        cpp_file.write(gradient_string)

    def write_grd_monomial(self, grd_file, index, monomial, variable_permutations):
        """
        Writes a single monomial to the 'poly-grd.maple' polynomial file.

        Args:
            cpp_file        - Local path to file to write in.
            index           - The index of this monomial.
            monomial        - The Monomial to write.
            variable_permutations - All permutations of the variables in this polynomial that can be used to permute
                    the Monomial.

        Returns:
            None.
        """
        grd_file.write("    p[{}] := ".format(index + 1))

        first_term = True
        for permutation in set(
                [permutation for permutation in monomial.permute(variable_permutations)]):

            if not first_term:
                grd_file.write("+")

            first_term = False

            first_factor = True

            for variable_index, degree in enumerate(permutation):

                if degree == 0:
                    continue

                for term in range(degree):

                    if not first_factor:
                        grd_file.write("*")

                    first_factor = False

                    grd_file.write("x{}".format(str(variable_index + 1).rjust(2, "0")))

        grd_file.write(":\n")

    def write_nogrd_monomial(self, nogrd_file, index, monomial, variable_permutations):
        """
        Writes a single monomial to the 'poly-nogrd.maple' polynomial file.

        Args:
            cpp_file        - Local path to file to write in.
            index           - The index of this monomial.
            monomial        - The Monomial to write.
            variable_permutations - All permutations of the variables in this polynomial that can be used to permute
                    the Monomial.

        Returns:
            None.
        """
        nogrd_file.write("    p[{}] := ".format(index + 1))

        first_term = True
        for permutation in set(
                [permutation for permutation in monomial.permute(variable_permutations)]):

            if not first_term:
                nogrd_file.write("+")

            first_term = False

            first_factor = True

            for variable_index, degree in enumerate(permutation):

                if degree == 0:
                    continue

                for term in range(degree):

                    if not first_factor:
                        nogrd_file.write("*")

                    first_factor = False

                    nogrd_file.write("x{}".format(str(variable_index + 1).rjust(2, "0")))

        nogrd_file.write(":\n")

    def write_cpp_closing(self, cpp_file, total_terms):
        """
        Writes the closing of the 'poly-direct.cpp' polynomial file.

        Args:
            cpp_file        - Local path to file to write in.
            total_terms     - Number of terms in the polynomial.

        Returns:
            None.
        """
        cpp_file.write("""    double energy(0);
    for(int i = 0; i < {}; ++i)
        energy += p[i]*a[i];

    return energy;

}}
}} // namespace mb_system""".format(total_terms))

    def write_cpp_grd_closing(self, cpp_file, total_terms):
        """
        Writes the closing of the 'poly-grd-direct.cpp' polynomial file.

        Args:
            cpp_file        - Local path to file to write in.
            total_terms     - Number of terms in the polynomial.

        Returns:
            None.
        """
        cpp_file.write("""    double energy(0);
    for(int i = 0; i < {}; ++i)
        energy += p[i]*a[i];

    return energy;

}}
}} // namespace mb_system""".format(total_terms))

    def write_grd_closing(self, grd_file, total_terms, number_of_variables):
        """
        Writes the closing of the 'poly-grd.maple' polynomial file.

        Args:
            cpp_file        - Local path to file to write in.
            total_terms     - Number of terms in the polynomial.
            number_of_variables - Number of variables in the polynomial.

        Returns:
            None.
        """
        arg_str = "["

        for index in range(number_of_variables):
            arg_str += "x" + str(index + 1).rjust(2, "0")
            if index != number_of_variables - 1:
                arg_str += ","
                if (index + 1) % 10 == 0:
                    arg_str += "\n         "
                else:
                    arg_str += " "

        arg_str += "]"
        grd_file.write("""
energy := 0;
for k from 1 by 1 to {} do
    energy := energy + a[k]*p[k]:
od:

args := {}:

energy := convert(energy, 'horner', args):

energy_proc := codegen[makeproc](energy, parameters = args):
codegen[cost](energy_proc);

xxx := codegen[GRADIENT](energy_proc, args, function_value = true):
codegen[cost](xxx);
xxx := codegen[optimize](xxx):
codegen[cost](xxx);

xxx := codegen[packargs](xxx, args, x):
xxx := codegen[optimize](xxx):

codegen[C](xxx, optimized, filename="poly-grd.c"):""".format(total_terms, arg_str))

    def write_nogrd_closing(self, nogrd_file, total_terms, number_of_variables):
        """
        Writes the closing of the 'poly-nogrd.maple' polynomial file.

        Args:
            cpp_file        - Local path to file to write in.
            total_terms     - Number of terms in the polynomial.
            number_of_variables - Number of variables in the polynomial.

        Returns:
            None.
        """
        arg_str = "["

        for index in range(number_of_variables):
            arg_str += "x" + str(index + 1).rjust(2, "0")
            if index != number_of_variables - 1:
                arg_str += ","
                if (index + 1) % 10 == 0:
                    arg_str += "\n         "
                else:
                    arg_str += " "

        arg_str += "]"
        nogrd_file.write("""
energy := 0;
for k from 1 by 1 to {} do
    energy := energy + a[k]*p[k]:
od:

args := {}:

energy := convert(energy, 'horner', args):

energy_proc := codegen[makeproc](energy, parameters = args):
codegen[cost](energy_proc);

xxx := codegen[optimize](energy_proc):
codegen[cost](xxx);

xxx := codegen[packargs](xxx, args, x):
xxx := codegen[optimize](xxx):

codegen[C](xxx, optimized, filename="poly-nogrd.c"):""".format(total_terms, arg_str))






