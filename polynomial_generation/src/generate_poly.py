import sys, os
import itertools

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../../")
import settings_reader
from exceptions import ParsingError, InvalidValueError, InconsistentValueError

def generate_poly(settings_file, input_file, order, output_path):
    
    # read the settings file
    settings = settings_reader.SettingsReader(settings_file)

    # open output files for polynomials, will automatically be closed by with open as syntax
    with open(output_path + "/poly.log", "w") as poly_log:
        
        # parse declaired fragments and variables from the input-file
        fragments, variables = parse_input(input_file)

        if len(fragments) == 1 and not settings.get("poly_generation", "accepted_terms") == "all":
            raise InconsistentValueError("number of fragments", "[poly_generation].accepted_terms", len(fragments), settings.get("poly_generation", "accepted_terms"), "when there is only 1 fragment, there are no intermolecular interactions, so you must use 'all' terms")

        # write number of fragments to log
        poly_log.write("<> fragments({}) <>\n".format(len(fragments)))
        poly_log.write("\n")

        # write each fragment to log
        for fragment in fragments:
            poly_log.write(fragment + "\n")
        poly_log.write("\n")

        # Parse the fragment names to name the atoms in the system

        # holds list of all atoms
        atom_names = []
        # number of atoms in each fragment
        number_per_fragment = []
        # index of first atom in atom_names of each fragment
        index_each_fragment = []
        # letter to be part of atom_name, signifies which fragment this atom is a part of starts at 'a', will be incremented to 'b', etc
        fragname = 97
        
        # loop over each fragment
        for fragment in fragments:

            # first check the formatting on this fragment
            if fragment[0].isdigit():
                raise InvalidValueError("fragment", fragment, "formatted such that the first character of every fragment is a letter")
            for i in range(1, len(fragment), 2):
                if not fragment[i].isdigit():
                    raise InvalidValueError("fragment", fragment, "formatted such that every letter is followed by a single-digit number, even if that number is 1")
            
            # counter to count how many atoms are in this fragment
            number_of_atoms = 0

            # loop thru all the atomic symbols in the fragment (skipping the numbers)
            for symbol in fragment[::2]:

                # loop once for each atom of this type as signified by the number immediately following the atomic symbol
                for k in range(0, int(fragment[fragment.index(symbol) + 1])):

                    # check if there are multiple atoms of this type in this fragment
                    if int(fragment[fragment.index(symbol) + 1]) > 1:

                        # add to atom_names in form Aa1
                        atom_names.append(symbol + str(k + 1) + chr(fragname))

                    else:

                        # add to atom_names in form Aa
                        atom_names.append(symbol + chr(fragname))
                    # increment number of atoms
                    number_of_atoms += 1

            # add index of first atom in this fragment to index_each_fragment
            index_each_fragment.append(sum(number_per_fragment))

            # add number of atoms in this fragment to number_per)fragment
            number_per_fragment.append(number_of_atoms)
            
            # increment fragname ('a' -> 'b', etc)
            fragname += 1

        # write number of atoms to log
        poly_log.write("<> atoms ({}) <>\n".format(len(atom_names)))
        poly_log.write("\n")

        # write each atom to log
        poly_log.write(":".join(atom_names) + "\n")
        poly_log.write("\n")

        # build the permutation group for each fragment

        fragment_permutations = []
        for frag_index, fragment in enumerate(fragments):
            permutations = list(make_permutations(fragment))
            
            for permutation in permutations:
                for index in range(len(permutation)):
                    permutation[index] += index_each_fragment[frag_index]

            fragment_permutations.append(permutations)

        # construct total permutations from the fragments, including permutations of like fragments

        atom_permutations = list(combine_permutations(fragments, fragment_permutations))

        # write permutation count to log file
        poly_log.write("<> permutations ({}) <>\n".format(len(atom_permutations)))
        poly_log.write("\n")

        # write each permutation to log file
        for permutation in atom_permutations:
            poly_log.write(":".join(str(x) for x in permutation) + "\n")
        poly_log.write("\n")

        # write variable count to log file
        poly_log.write("<> variables ({}) <>\n".format(len(variables))) 
        poly_log.write("\n")

        # generate the variable permutations
        variable_permutations = list(make_variable_permutations(variables, atom_permutations, atom_names))

        # make the .cpp vars file
        with open(output_path + "/vars.cpp", "w") as vars_file:
            write_variable_file(vars_file, variables);

        # write each variable to log file
        for index, variable in zip(range(len(variables)), variables):
            poly_log.write("{:>2} : {:>3}({:>1}) <===> {:>3}({:>1}) : {}\n".format(index, variable.atom1_name, variable.atom1_fragment, variable.atom2_name, variable.atom2_fragment, variable.category))
        poly_log.write("\n")

        # generate the monomials

        total_terms = 0

        # this list will be filled so that the first item is all the 1st degree monomials, second item is all the 2nd degree monomials, etc
        total_monomials = []

        # loop thru every degree in this polynomial
        for degree in range(1, order + 1):

            # header for this degree
            poly_log.write("<> {} degree <>\n".format(degree))
            poly_log.write("\n")

            # get all the monomials of the current degree
            monomials = list(generate_monomials(len(variables), degree))

            # log number of possible monomials
            poly_log.write("{} possible {} degree monomials\n".format(len(monomials), degree))

            if settings.get("poly_generation", "accepted_terms") == "all":
                accepted_monomials = monomials
            elif settings.get("poly_generation", "accepted_terms") == "partly-inter":
                accepted_monomials = list(eliminate_purely_intramolecular_monomials(monomials, variables))
            elif settings.get("poly_generation", "accepted_terms") == "purely-inter":
                accepted_monomials = list(eliminate_partly_intramolecular_monomials(monomials, variables))
            else:
                raise InvalidValueError("[poly_generation][accepted_terms]", settings.get("poly_generation", "accepted_terms"), "one of 'all', 'partly-inter', or 'inter'")

            # filter out redundant monomials (that are a permutation of eachother)
            accepted_monomials = list(eliminate_redundant_monomials(accepted_monomials, variable_permutations))

            # log number of accpeted terms
            poly_log.write("{} <<== accepted {} degree terms\n".format(len(accepted_monomials), degree))

            poly_log.write("\n")

            # update the total number of terms
            total_terms += len(accepted_monomials)

            # add the monomials of the current degree to the list of all monomials
            total_monomials.append(accepted_monomials)

        # log the total number of terms
        poly_log.write(" Total number of terms: {}\n".format(total_terms))

        # write the header file
        with open(output_path + "/poly-model.h", "w") as header_file:
            write_header_file(header_file, total_terms, len(variables))

        # open the three files we will now write to
        with open(output_path + "/poly-direct.cpp", "w") as cpp_file, open(output_path + "/poly-grd.maple", "w") as grd_file, open(output_path + "/poly-nogrd.maple", "w") as nogrd_file:
            # write the opening for the cpp file
            write_cpp_opening(cpp_file, total_terms, len(variables))

            # keeps track of what index in a list of all monomials the current monomial would occupy
            monomial_index = 0

            # loop thru every degree in this polynomial
            for degree in range(1, order + 1):
                for monomial in total_monomials[degree - 1]:
                    write_cpp_monomial(cpp_file, monomial_index, monomial, variable_permutations)
                    write_grd_monomial(grd_file, monomial_index, monomial, variable_permutations)
                    write_nogrd_monomial(nogrd_file, monomial_index, monomial, variable_permutations)
                    monomial_index += 1

                cpp_file.write("\n")
                grd_file.write("\n")

            write_cpp_closing(cpp_file, total_terms)
            write_grd_closing(grd_file, total_terms, len(variables))
            write_nogrd_closing(nogrd_file, total_terms, len(variables))


def parse_input(input_path):
    """
    Reads the add_molecule, and add_fragments statements from the input file
    """

    fragments = []
    variables = []

    with open(input_path, "r") as input_file:

        # loop thru all lines in the input file
        for line in input_file:

            # check if this line is an add_molecule statement
            if line[:12] == "add_molecule":

                # parse line into a fragment
                fragments.append(line[line.index("['") + 2:line.index("']")])

            # check if this line is an add_variable statement
            elif line[:12] == "add_variable":

                # parse line into a variable
                variables.append(Variable(line))

            # if line is neither a add_molecule statement, add_variable statement, or a blank line, it is invalid
            elif line != "\n":
                raise ParsingError(input_path, "Unrecongnized line format: '{}'".format(line))

    return fragments, variables

class Variable(object):
    """
    Holds all information relevent to a single variable
    """
    def __init__(self, line):
        # remove all the extra information before and after the brackets
        line = line[line.index("[") + 1:line.index("]")]

        # parse the line into this Variable's fields
        self.atom1_name, self.atom1_fragment, self.atom2_name, self.atom2_fragment, self.category = line.replace("'", "").replace(" ", "").split(",")
    

def make_permutations(fragment):
    """
    Takes in a fragment and constructs all permutations of that fragment
    """

    try:
        # read the 2nd letter of the fragment, the atom count. (For instance 2 in A2B3)
        count = int(fragment[1])
    except IndexError:
        # if the fragment string is empty, then this fragment has no permutations
        yield []
        return

    # get all possible permutations of the number of atoms indicated by count
    permutations = itertools.permutations(range(count))
    
    # loop thru each of the permutations
    for permutation in permutations:
        
        # yield from the generator created by concatinating the permutation of the first atom type in the fragment with the permutations created by a recursive call on the rest of the fragment
        # increase each item in the lists given by the recursive call equal to the number of atoms of the first type in the fragment
        yield from (list(permutation) + [count + i for i in perm] for perm in make_permutations(fragment[2:]))

def combine_permutations(fragments, fragment_permutations):
    """
    Takes a list of fragments ("A1B2", "A1B1C1", etc) and a list of the permutations of those fragments and yields permutations of the molecules as a whole
    """
    
    # if there are no fragments, then there are no permutations for this molecule
    if len(fragments) == 0:
        yield []
        return

    # loop thru every permutation of the first fragment
    for permutation in fragment_permutations[0]:

        # yield from the first fragment's permutation concatinated with the permutations of the rest of the molecule
        yield from (permutation + perm for perm in combine_permutations(fragments[1:], fragment_permutations[1:]))

    # loop thru every other fragmanet in the molecule
    for i in range(1, len(fragments)):
        
        # check if the current fragment (0) is the same as another one
        if fragments[i] == fragments[0]:

            # loop thru every permutation of the fragment that is the same as the current one
            for permutation in fragment_permutations[i]:

                # yield from that fragment's permutation concatinated with the permuations of the rest of the molecule with the current fragment (0) swapped into the position of the fragment it is the same as.
                yield from (permutation + perm for perm in combine_permutations(fragments[1:i] + [fragments[0]] + fragments[i + 1:], fragment_permutations[1:i] + [fragment_permutations[0]] + fragment_permutations[i + 1:]))

def make_variable_permutations(variables, atom_permutations, atom_names):
    """
    Constructs the variable permutations from the atom permutations
    """

    # there will be 1 variable permutation for every atom permutation
    for atom_permutation in atom_permutations:

        # initialize the new variable permutation
        variable_permutation = []

        # each variable permutation has length = len(variables) where each value is a value of a variable which should be switched with this index to create this permutation
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
                if new_atom1 == new_variable.atom1_name + new_variable.atom1_fragment and new_atom2 == new_variable.atom2_name + new_variable.atom2_fragment:
                    new_index = new_variable_index
                    break
                if new_atom2 == new_variable.atom1_name + new_variable.atom1_fragment and new_atom1 == new_variable.atom2_name + new_variable.atom2_fragment:
                    new_index = new_variable_index
                    break

            # a variable should always be found, if one is not found new_index being -1 signifies a problem
            if new_index == -1:
                print("Something went wrong :(")

            # append the index of the permutated variable to variable permutation
            variable_permutation.append(new_index)

        yield variable_permutation

            

def generate_monomials(number_of_vars, degree):
    """
    Given a number of variables and a degree, generates all possible monomials of that degree.

    yeilds arrays of length number_of_vars with entries adding to degree
    """

    # if number of vars is 1, then the only possible monomial is a monomial with 1 variable with the given degree
    if number_of_vars == 1:
        yield [degree]
        return

    # if degree is 0, then the only possible monomial is a monomial with all 0 degrees with the given number of variables
    if degree == 0:
        yield [0 for i in range(number_of_vars)]
        return

    # loop over all possible first non-zero terms
    for v in range(number_of_vars):

        # loop over all possible degree values of the first non-zero term
        for d in range(1, degree + 1):

            # yield from the list created by all zeros before the first non-zero term, then the first non-zero term, then each result of the recursive call on all terms after the first non-zero term with degree equal to degree minus the degree of the first non-zero term
            yield from ([0 for i in range(v)] + [d] + monomial for monomial in generate_monomials(number_of_vars - v - 1, degree - d))

def eliminate_redundant_monomials(monomials, variable_permutations):
    accepted_monomials = monomials[:]
    index = 0
    while(True):
        try:
            monomial = accepted_monomials[index]
        except IndexError:
            break
        permutated_monomials = set([tuple(permutation) for permutation in permute_monomial(monomial, variable_permutations)])
        for permutated_monomial in permutated_monomials:
            try:
                accepted_monomials.remove(list(permutated_monomial))
            except ValueError:
                pass

        accepted_monomials.insert(0, monomial)
        index += 1

    yield from (x for x in accepted_monomials)
    
def permute_monomial(monomial1, variable_permutations):
    for variable_permutation in variable_permutations:

        monomial1_permutation = [0 for i in monomial1]

        for index, degree in zip(range(len(monomial1)), monomial1):

            if degree == 0:
                continue
            
            new_index = variable_permutation.index(index)
            monomial1_permutation[new_index] = monomial1[index]

        yield monomial1_permutation
            
def eliminate_purely_intramolecular_monomials(monomials, variables):
    for monomial in monomials:
        if not is_purely_intramolecular(monomial, variables):
            yield monomial

def eliminate_partly_intramolecular_monomials(monomials, variables):
    for monomial in monomials:
        if not is_partly_intramolecular(monomial, variables):
            yield monomial

def is_purely_intramolecular(monomial, variables):
    for index, degree in enumerate(monomial):
        if degree > 0 and not variables[index].category[:7] == "x-intra":
            return False
    return True

def is_partly_intramolecular(monomial, variables):
    for index, degree in enumerate(monomial):
        if degree > 0 and variables[index].category[:7] == "x-intra":
            return True
    return False

def write_variable_file(variable_file, variables):
    # variable header comment
    variable_file.write("\n// <-> variables <->\n\n")
    
    # loop thru each variable
    for index, variable in enumerate(variables):
        
        # write the variable to vars.cpp
        variable_file.write("    x[{}] = @VAR@-|{}|({}{}{}{};\n".format(index, variable.category, variable.atom1_name, variable.atom1_fragment, variable.atom2_name, variable.atom2_fragment))

def write_header_file(header_file, total_terms, number_of_variables):
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

public:
    unsigned report_nvars(){{ return n_vars; }};
    unsigned report_size(){{ return size; }};
}};

}} // namespace mb_system

#endif // POLY_MODEL_H""".format(total_terms, number_of_variables))

def write_cpp_opening(cpp_file, total_terms, number_of_variables):
    cpp_file.write("""#include "poly-model.h"

namespace mb_system {{

double poly_model::eval_direct(const double a[{0}], const double x[{1}])
{{
    double p[{0}];
""".format(total_terms, number_of_variables))

def write_cpp_monomial(cpp_file, index, monomial, variable_permutations):
    cpp_file.write("    p[{}] = ".format(index))

    first_term = True
    for permutation in set([tuple(permutation) for permutation in permute_monomial(monomial, variable_permutations)]):

        if not first_term:
            cpp_file.write(" + ")

        first_term = False;

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

def write_grd_monomial(grd_file, index, monomial, variable_permutations):
    grd_file.write("    p[{}] := ".format(index))

    first_term = True
    for permutation in set([tuple(permutation) for permutation in permute_monomial(monomial, variable_permutations)]):

        if not first_term:
            grd_file.write("+")

        first_term = False;

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

def write_nogrd_monomial(nogrd_file, index, monomial, variable_permutations):
    nogrd_file.write("    p[{}] := ".format(index))

    first_term = True
    for permutation in set([tuple(permutation) for permutation in permute_monomial(monomial, variable_permutations)]):

        if not first_term:
            nogrd_file.write("+")

        first_term = False;

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

def write_cpp_closing(cpp_file, total_terms):
    cpp_file.write("""    double energy(0);
    for(int i = 0; i < {}; ++i)
        energy += p[i]*a[i];

    return energy;

}}
}} // namespace mb_system""".format(total_terms))

def write_grd_closing(grd_file, total_terms, number_of_variables):
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

def write_nogrd_closing(nogrd_file, total_terms, number_of_variables):
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

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python generate_poly.py settings.ini input_file order output_path")
        sys.exit(1)
    generate_poly(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
