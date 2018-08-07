import sys
import configparser
import itertools

def generate_poly(settings, input_file, order, output_path):
    
    # create configparser
    config = configparser.SafeConfigParser(allow_no_value=False)
    config.read(settings)

    # open output files for polynomials, will automatically be closed by with open as syntax
    with open(output_path + "/vars.cpp", "w") as out_vars, open(output_path + "/poly-direct.cpp", "w") as out_cpp, open(output_path + "/poly-nogrd.maple", "w") as out_maple_nogrd, open(output_path + "/poly-grd.maple", "w") as out_maple_grd, open(output_path + "/poly.log", "w") as poly_log, open(input_file, "r") as in_poly:
        
        # read the add_molecule definitions to get a list of the fragments
        fragments = parse_fragments(in_poly)

        # write number of fragments to log
        poly_log.write("<> fragments({}) <>\n".format(len(fragments)))
        poly_log.write("\n")

        # write each fragment to log
        for fragment in fragments:
            poly_log.write(fragment + "\n")
        poly_log.write("\n")

        # Parse the molecule names to determine and name the atoms in the system

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

            # counter to count how many atoms are in this fragment
            number_of_atoms = 0

            # loop thru all the atomic symbols in the fragment (skipping the numbers)
            for symbol in fragment[::2]:
                # loop once for each atom of this type as signified by the number immediately following the atomic symbol
                for k in range(0, int(fragment[fragment.index(symbol) + 1])):
                    # check if there are multiple atoms of this type in this fragment
                    if int(fragment[fragment.index(symbol) + 1]) > 1:
                        # add to atom_names in form Aa1
                        atom_names.append(symbol + chr(fragname) + str(k + 1))
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

        # calculate total number of atoms in all fragments
        num_atoms = len(atom_names)

        # write number of atoms to log
        poly_log.write("<> atoms ({}) <>\n".format(num_atoms))
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

        # read variables from input file

        variables = list(parse_variables(in_poly))

        # write variable count to log file
        poly_log.write("<> variables ({}) <>\n".format(len(variables))) 
        poly_log.write("\n")

        # write each variable to log file
        for index, variable in zip(range(len(variables)), variables):
            poly_log.write("{:>2} : {:>3}({:>1}) <===> {:>3}({:>1}) : {}\n".format(index, variable.atom1_name, variable.atom1_fragment, variable.atom2_name, variable.atom2_fragment, variable.category))
        poly_log.write("\n")

        # generate the monomials

        for degree in range(1, order + 1):
            poly_log.write("<> {} degree <>\n".format(degree))
            poly_log.write("\n")
            monomials = list(generate_monomials(len(variables), degree))
            poly_log.write("{} possible {} degree monomials\n".format(len(monomials), degree))
        
        
def parse_fragments(input_file):
    """
    Reads the add_molecule statements from the input file
    """

    fragments = []

    # get a line from the input file
    line = input_file.readline()

    # continue reading lines as long as they begin with add_m
    while line[:5] == "add_m":
        
        # parse the fragment from the line
        fragments.append(line[line.index("['") + 2:line.index("']")])

        # get a new line from the input file
        line = input_file.readline()

    return fragments

def parse_variables(input_file):
    """
    Reads the add variable statements from the input file
    """

    # loop thru all the rest of the lines in the input file, constructing variable objects
    for line in input_file.readlines():
        yield Variable(line)
        
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

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python generate_poly.py settings.ini input_file order output_path")
        sys.exit(1)
    generate_poly(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
