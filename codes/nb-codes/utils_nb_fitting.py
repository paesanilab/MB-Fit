from potential_fitting.polynomials import FragmentParser
from potential_fitting.exceptions import InconsistentValueError

def get_atom_types(fragment):
    """
    Returns a list with the atom types and the number of atoms of that type. For A1B3, returns ["A",1,"B",3]

    Args:
        fragment              - The fragment you want to decompose
    """

    fragment_parser = FragmentParser(fragment, 'a')

    atom_list = []

    for atom_type in fragment_parser.get_atom_and_virtual_site_types():
        atom_list.append(atom_type.get_type())
        atom_list.append(atom_type.get_count())

    return atom_list

def get_nonbonded_pairs(vsites, mon_types_a, mon_types_b = None):
    """
    Returns a list with all the pairs belonging to real atoms in the order they appear in both monomers, without repeating them. ['A', 1, 'B', 4, 'C', 2] will return ['AA', 'AB', 'AC', 'BB', 'BC', 'CC'], while if two arguments for the types are passed it will do all the pairs.

    Args:
        vsites                - Virtual site labels
        mon_types_a           - List with the atom symmetry and the number of atoms of that symmetry, e.g. [[A,1,B,1],[C,2]] for A1B2_C2 for the first monomer
        mon_types_b           - List with the atom symmetry and the number of atoms of that symmetry, e.g. [[A,1,B,1],[C,2]] for A1B2_C2 for the second monomer
    """

    pairs = []
    if mon_types_b is None:
        for i in range(0,len(mon_types_a),2):
            for j in range(i,len(mon_types_a),2):
                this_pair = "".join(sorted([mon_types_a[i],mon_types_a[j]]))
                if not this_pair in pairs and not any(x in this_pair for x in vsites):
                    pairs.append(this_pair)
    else:
        for i in range(0,len(mon_types_a),2):
            for j in range(0,len(mon_types_b),2):
                this_pair = "".join(sorted([mon_types_a[i],mon_types_b[j]]))
                if not this_pair in pairs and not any(x in this_pair for x in vsites):
                    pairs.append(this_pair)
    return pairs
    
    

def read_poly_in(poly_in, vsites, var_intra, var_inter, var_virtual_sites):
    """
    Reads the poly.in file, and returns a list of variables, the intramolecular pairs and the intermolecular pairs.

    Args:
        poly_in               - The polynomial input file
        vsites                - Virtual site labels
        var_intra             - The intramolecular variable type
        var_inter             - The intermolecular variable type
        var_virtual_sites     - The virtual sites variable type
    """

    variables = []
    intra_poly_pairs = []
    inter_poly_pairs = []
    with open(poly_in, "r") as input_file:
        for line in input_file:
            if line.startswith("add_variable"):
                args = line[line.index('[') + 1:line.index(']')].replace("'", "").replace(" ", "").split(",")
    
                variables.append(args)
    
                pair = "".join(sorted([args[0][0], args[2][0]]))
                
                # Check if there is a virtual site involved
                has_vsites = any(x in pair for x in vsites)
    
                if args[4].startswith("x-intra"):
                    if (args[1] == args[3]):
                        if has_vsites:
                            variables[-1].append(var_virtual_sites)
                        else:
                            variables[-1].append(var_intra)
                        if pair not in intra_poly_pairs:
                            intra_poly_pairs.append(pair)
                    else:
                        # variable is labeled as intra but the two fragment labels are different. UH OH!
                        raise InconsistentValueError('fragment id of first atom', 'fragment id of second atom', args[1], args[3], 'x-intra variables must have the same fragment index; something is wrong in {}'.format(poly_in))
                else:

                    if has_vsites:
                        variables[-1].append(var_virtual_sites)
                    else:
                        variables[-1].append(var_inter)
                    if pair not in inter_poly_pairs:
                        inter_poly_pairs.append(pair)
    return variables, intra_poly_pairs, inter_poly_pairs

def get_non_linear_parameters(variables):
    """
    From the variables it returns all the non linear parameters, that will be used in the polynomial evaluation

    Args:
        variables              - List of lists with the information in the poly.in file
    """

    intra_nl_params = []
    inter_nl_params = []
    nl_params_ordered = []
    for i in range(len(variables)):
        # Obtain pair from variables
        pair = "".join(sorted([get_atom_types(variables[i][0])[0], get_atom_types(variables[i][2])[0]]))
        
        # Check if pair functional form has d0. If so, add it.
        if "0" in variables[i][-1]:
            nlp_to_use = ['k','d']
        else:
            nlp_to_use = ['k']

        # See if pair is intra or inter
        if variables[i][4].startswith("x-intra"):
            # loop over different constants
            for nl_p in nlp_to_use:
                nl_constant = "{}_intra_{}".format(nl_p, pair)
                if not nl_constant in intra_nl_params:
                    intra_nl_params.append(nl_constant)
                    nl_params_ordered.append(nl_constant)
        else:
            for nl_p in nlp_to_use:
                nl_constant = "{}_{}".format(nl_p, pair)
                if not nl_constant in inter_nl_params:
                    inter_nl_params.append(nl_constant)
                    nl_params_ordered.append(nl_constant)
            
    return intra_nl_params,inter_nl_params, nl_params_ordered    

def get_list_of_numeric_pairs(prefix,number_of_monomers):
    """
    Helper function that returns the combinations of pairs within the number of monomers prepended with a string. This is to get the switches and the distances; for n monomers (if prefix is d) [d12, d13,...,d1n,d23, d24...]

    Args:
        number_of_monomers     - Number of monomers in the system
    """

    my_labels = []
    for i in range(1,number_of_monomers):
        for j in range(i+1,number_of_monomers + 1):
            my_labels.append([prefix + str(i) + str(j), i, j])

    return my_labels








