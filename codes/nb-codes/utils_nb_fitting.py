def get_atom_types(fragment):
    atom_list = []
    was_digit = False
    current_text = ""
    for character in fragment:
        if not character.isdigit():
            if was_digit:
                atom_list.append(int(current_text))
                current_text = character
                was_digit = False
            else:
                current_text += character
        else:
            if was_digit:
                current_text += character
            else:
                atom_list.append(current_text)
                was_digit = True
                current_text = character
            
    # At this point, only the last number is missing
    atom_list.append(int(current_text))
    return atom_list

def get_nonbonded_pairs(vsites, mon_types_a, mon_types_b = None):
    # Case in which we only pass one monomer (for 1b)
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
                        raise ValueError # figure out what exception to raise
                else:
                    if has_vsites:
                        variables[-1].append(var_virtual_sites)
                    else:
                        variables[-1].append(var_inter)
                    if pair not in inter_poly_pairs:
                        inter_poly_pairs.append(pair)
    return variables, intra_poly_pairs, inter_poly_pairs

def get_non_linear_parameters(variables):
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
    # returns, for n monomers (if prefix is d) [d12, d13,...,d1n,d23, d24...]
    my_labels = []
    for i in range(1,number_of_monomers):
        for j in range(i+1,number_of_monomers + 1):
            my_labels.append([prefix + str(i) + str(j), i, j])

    return my_labels








