from potential_fitting.utils import SettingsReader
from potential_fitting.exceptions import InvalidValueError, InconsistentValueError

def generate_input_poly(settings_file, poly_input):
    """
    Generates the input for the polynomial generator

    Will be file with name like A1B2C2.in, etc as specified in the settings file
    """

    # create SettingsReader currently unused, but in future could be used for stuff like adding filters
    settings = SettingsReader(settings_file)

    # get the string representation of the molecule by getting the A1B2C2, etc part of the name
    try:
        molecule = os.path.splitext(poly_input)[0].split("/")[-1]
    except IndexError:
        raise InvalidValueError("poly_input", poly_input, "in format somepath/A1B2_A1B2.in")
    # file will be automatically closed after block by with open as ... syntax
    with open(poly_input, "w") as poly_in_file:

        # split along _ to seperate fragments
        fragments = molecule.split('_')

        # loop thru each fragment and add an add_molecule line at top of file
        for fragment in fragments:
            poly_in_file.write("add_molecule['" + fragment + "']\n")

        # these letters mark virtual sites instead of atoms
        virtual_sites = ['X', 'Y', 'Z']

        # loop thru each type list and generate each fragment's set of uniquely named atoms
        sets_list = []

        # starts at 97 ('a') for first fragment, incremeneted for each fragment
        letter = 97

        # loop thru every fragment
        for fragment in fragments:

            # first check the formatting on this fragment
            if fragment[0].isdigit():
                raise InvalidValueError("poly_input", poly_input, "formatted such that the first character of every fragment is a letter")
            for i in range(1, len(fragment), 2):
                if not fragment[i].isdigit():
                    raise InvalidValueError("poly_input", poly_input, "formatted such that every letter is followed by a single-digit number, even if that number is 1")
            frag_set = []
            
            # loop thru every other character in fragment, get the letters (A, B, C) but not numbers (1, 2, 3)
            for i in range(0, len(fragment), 2):

                # check if there is only one of this atom in this fragment by looking at the character after the letter, which is a number signifying how many of that atom are in this fragment
                if (int(fragment[i + 1]) == 1):

                    # check if this atom is not a virtual site
                    if not fragment[i] in virtual_sites:

                        # add it to the fragment's set without a number in format A__a
                        frag_set.append(fragment[i] + '_' +  '_' + chr(letter))

                # otherwise there is more than one of this atom in this fragment
                else:

                    # counter to count how many of this atom are in the fragment
                    num_atoms = 1

                    # loop thru as many atoms as there are of this type in this fragment
                    for j in range(int(fragment[i + 1])):
                        
                        # check if this atom is not a virtual site
                        if not fragment[i] in virtual_sites:

                            # add it to the fragment's set with a number in format A_1_a
                            frag_set.append(fragment[i] + '_' + str(num_atoms) + '_' + chr(letter))
                            num_atoms += 1

            # loop thru every character in fragment again, this time to look for virtual sites
            for i in range(0, len(fragment), 2):

                # counter to count how many of this virtual site are in the fragment
                num_atoms = 1

                # check if atom is a virtual site
                if fragment[i] in virtual_sites:

                    # look one character ahead to see how many of these virtual sites are in the fragment
                    for j in range(int(fragment[i + 1])):

                        # add virtual site to set in format X_1_a
                        frag_set.append(fragment[i] + '_' + str(num_atoms) + '_' + chr(letter))
                        num_atoms += 1

            # add this fragment's set to list of sets
            sets_list.append(frag_set)

            # increment the letter, so next fragment has next lowercase letter alphabetically
            letter += 1

        # write new line to file
        poly_in_file.write('\n')
        
        # loop thru each fragment and write its variables to the file
        for frag_set in sets_list:

            # loop thru every combination of atoms in this fragment
            for i in range(0, len(frag_set) - 1):
                for j in range(i + 1, len(frag_set)):

                    # split the first atom along _ to seperate it: "A_1_a" -> ["A", "1", "a"] or "A__a" -> ["A", "", "a"]
                    split_atom1 = frag_set[i].split('_')
                    # split the second atom along _ to seperate it: "A_1_a" -> ["A", "1", "a"] or "A__a" -> ["A", "", "a"]
                    split_atom2 = frag_set[j].split('_')

                    # combine first element of each split atom in alphabetical order to get something like "AA", or "AB"
                    combined_letters = ''.join(sorted(split_atom1[0] + split_atom2[0]))

                    # check if neither atom is a virtual site
                    if not split_atom1[0] in virtual_sites and not split_atom2[0] in virtual_sites:

                        # write valiable to file in proper format
                        poly_in_file.write("add_variable['" + split_atom1[0] + split_atom1[1] + "', '" + split_atom1[2] + "', '" + split_atom2[0] + split_atom2[1] + "', '" + split_atom2[2] + "', 'x-intra-" + combined_letters + "']\n")
        
        # loop thru every pair of fragments and add write variables to file for interfragmental interactions
        for frag_set1 in sets_list:
            # loop thru all frag_set2 that come after frag_set1
            for frag_set2 in sets_list[sets_list.index(frag_set1) + 1:]:

                # loop thru every comination of atoms from the frag_sets
                for i in range(0,len(frag_set1)):
                    for j in range(0,len(frag_set2)):
                        
                        # split the first atom along _ to seperate it: "A_1_a" -> ["A", "1", "a"] or "A__a" -> ["A", "", "a"]
                        split_atom1 = frag_set1[i].split('_')
                        # split the second atom along _ to seperate it: "A_1_a" -> ["A", "1", "a"] or "A__a" -> ["A", "", "a"]
                        split_atom2 = frag_set2[j].split('_')
                        
                        # combine first element of each split atom in alphabetical order to get something like "AA", or "AB"
                        combined_letters = ''.join(sorted(split_atom1[0] + split_atom2[0]))
                        poly_in_file.write("add_variable['" + split_atom1[0] + split_atom1[1] + "', '" + split_atom1[2] + "', '" + split_atom2[0] + split_atom2[1] + "', '" + split_atom2[2] + "', 'x-" + combined_letters + "']\n")

        # newline between variables and filters
        poly_in_file.write("\n")

        # add filter based on the filtering setting in settings.ini
        polynomial_filtering = settings.get("poly_generation", "accepted_terms", "all" if len(fragments) == 1 else "purely-inter")

        # make sure the user hasn't chosed intermolecular terms only with a monomer
        if len(fragments) == 1 and not polynomial_filtering == "all":
            raise InconsistentValueError("number of fragments", "[poly_generation].accepted_terms", len(fragments), polynomial_filtering, "when there is only 1 fragment, there are no intermolecular interactions, so you must use 'all' terms")

        if polynomial_filtering == "purely-inter":
            # this filter filters out all terms that have any intra-molecular components
            poly_in_file.write("add_filter['degree', 'x-intra-**', '1+', '*']")
        elif polynomial_filtering == "partly-inter":
            # this filter filters out all terms that have no inter-molecular components
            poly_in_file.write("add_filter['not', 'degree', 'x-**', '1+', '*']")

        elif polynomial_filtering != "all":
            raise InvalidValueError("[poly_generation][accepted_terms]", polynomial_filtering, "one of 'all', 'partly-inter', or 'purely-inter'")
