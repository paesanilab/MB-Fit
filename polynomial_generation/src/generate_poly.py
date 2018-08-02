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
        for fragment in fragments:
            permutations = list(make_permutations([], fragment, 0))
            
            for permutation in permutations:
                for index in range(len(permutation)):
                    permutation[index] += index_each_fragment[fragments.index(fragment)]

            fragment_permutations.append(permutations)

        # construct total permutations from the fragments, including permutations of like fragments

        atom_permutations = list(combine_permutations([], fragments, fragment_permutations))

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
            monomials = list(generate_monomials(list(0 for x in range(len(variables))), degree, 0))
            poly_log.write("{} possible {} degree monomials\n".format(len(monomials), degree))
        
        

def parse_fragments(input_file):
    fragments = []
    line = input_file.readline()
    while line[:5] == "add_m":
        fragments.append(line[line.index("['") + 2:line.index("']")])
        line = input_file.readline()

    return fragments

def parse_variables(input_file):
    for line in input_file.readlines():
        yield Variable(line)
        
class Variable(object):
    def __init__(self, line):
        line = line[line.index("[") + 1:line.index("]")]
        self.atom1_name, self.atom1_fragment, self.atom2_name, self.atom2_fragment, self.category = line.replace("'", "").replace(" ", "").split(",")
    

def make_permutations(existing_permutation, fragment, index):
    # 2nd letter of fragment is next atom count
    try:
        count = int(fragment[1])
    except IndexError:
        yield existing_permutation
        raise StopIteration

    permutations = itertools.permutations(range(count))
    
    for permutation in permutations:
        yield from make_permutations(existing_permutation + [i + index for i in permutation], fragment[2:], index + count)

def combine_permutations(existing_permutation, fragments, fragment_permutations):

    if len(fragments) == 0:
        yield existing_permutation
        raise StopIteration

    for permutation in fragment_permutations[0]:
        yield from combine_permutations(existing_permutation + permutation, fragments[1:], fragment_permutations[1:])

    for i in range(1, len(fragments)):
        if fragments[i] == fragments[0]:
            for permutation in fragment_permutations[i]:
                yield from combine_permutations(existing_permutation + permutation, fragments[1:i] + [fragments[0]] + fragments[i + 1:], fragment_permutations[1:i] + [fragments[0]] + fragments[i + 1:])

def generate_monomials(existing_powers, degree, min_index):
    if degree == 0:
        yield existing_powers
        raise StopIteration

    for i in range(min_index, len(existing_powers)):
        new_powers = existing_powers[:]
        new_powers[i] += 1
        yield from generate_monomials(new_powers, degree - 1, i)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python generate_poly.py settings.ini input_file order output_path")
        sys.exit(1)
    generate_poly(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
