import sys
import os
import configparser

def generate_input_poly(settings):
    """
    Generates the input for the polynomial generator

    Will be file with name like A1B2C2.in, etc
    """

    # create configparser
    config = configparser.SafeConfigParser(allow_no_value=False)
    config.read(settings)

    # get the string representation of the molecule by getting the A1B2C2, etc part of the name
    molecule = os.path.splitext(config["files"]["poly_in_path"])[-2]

    with open(config["files"]["poly_in_path"], "w") as poly_in_file:
        # split along _ to seperate fragments
        number_of_fragments = len(molecule.split('_'))

        # get list of fragments
        fragments = molecule.split('_')

        # loop thru each fragment and add an add_molecule line at top of file
        for fragment in fragments:
            poly_in_file.write("add_molecule['" + fragment + "']\n")

        # these letters mark virtual sites instead of atoms
        virtual_sites = ['X', 'Y', 'Z']

        # loop thru each fragment and make a types list out of it
        types_list = []
        for fragment in fragments:

            # create types list
            types = list(fragment)
            # add types list to list of types lists
            types_list.append(types)

        # loop thru each type list and generate each fragment's "set"
        sets_list = []
        letter = 97
        for types in types_list:
            frag_set = []
            for i in range(0, len(types), 2):
                n = 1
                if (int(types[i + 1]) == 1):
                    if not types[i] in virtual_sites:
                        frag_set.append(types[i] + '_' +  '_' + chr(letter))
                else:
                    for j in range(int(types[i + 1])):
                        if not types[i] in virtual_sites:
                            frag_set.append(types[i] + '_' + str(n) + '_' + chr(letter))
                            n = n + 1

            for i in range(0, len(types), 2):
                n = 1
                for j in range(int(types[i + 1])):
                    if types[i] in virtual_sites:
                        frag_set.append(types[i] + '_' + str(n) + '_' + chr(letter))
                        n = n + 1
            # add this fragment's set to list of sets
            sets_list.append(frag_set)
            letter += 1
    
        # write new line to file
        poly_in_file.write('\n')
        
        # loop thru each fragment and write its variables to the file
        for frag_set in sets_list:
            for i in range(0, len(frag_set) - 1):
                for j in range(i + 1, len(frag_set)):
                    ti = frag_set[i].split('_')
                    tj = frag_set[j].split('_')
                    t = ''.join(sorted(ti[0] + tj[0]))
                    if not ti[0] in virtual_sites and not tj[0] in virtual_sites:
                        poly_in_file.write("add_variable['" + ti[0] + ti[1] + "', '" + ti[2] + "', '" + tj[0] + tj[1] + "', '" + tj[2] + "', 'x-intra-" + t + "']\n")
        
        # loop thru every pair of fragments and add variables for interfragmental interactions
        for frag_set1 in sets_list:
            for frag_set2 in sets_list[sets_list.index(frag_set1) + 1:]:
                
                for i in range(0,len(frag_set1)):
                    for j in range(0,len(frag_set2)):
                        ti = frag_set1[i].split('_')
                        tj = frag_set2[j].split('_')
                        t = ''.join(sorted(ti[0] + tj[0]))
                        poly_in_file.write("add_variable['" + ti[0] + ti[1] + "', '" + ti[2] + "', '" + tj[0] + tj[1] + "', '" + tj[2] + "', 'x-" + t + "']\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python generate_input_poly.py settings.ini")
        sys.exit(1)
    generate_input_poly(sys.argv[1])
