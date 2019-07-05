# external package imports
import itertools

# absolute module imports
from potential_fitting.utils import constants

class MoleculeInParser(object):

    def __init__(self, molecule_in):
        self.fragment_parsers = []
        print("molecule_init", molecule_in)
        frag_id = 'a'

        for fragment_in in molecule_in.split("_"):
            print("looping...")
            self.fragment_parsers.append(FragmentParser(fragment_in, frag_id))

            frag_id = chr(ord(frag_id) + 1)

    def get_molecule_in(self):
        return "_".join([fragment_parser.get_fragment_in() for fragment_parser in self.fragment_parsers])

    def get_intra_molecular_variables(self):
        for fragment_parser in self.get_fragments():

            yield from fragment_parser.get_intra_molecular_variables()

    def get_inter_molecular_variables(self):
        for fragment_parser_index, fragment_parser1 in enumerate(self.get_fragments()):

            for fragment_parser2 in list(self.get_fragments())[fragment_parser_index + 1:]:

                yield from fragment_parser1.get_inter_molecular_variables(fragment_parser2)

    def get_variables(self):
        yield from self.get_intra_molecular_variables()
        yield from self.get_inter_molecular_variables()

    def get_fragments(self):

        for fragment_parser in self.fragment_parsers:

            yield fragment_parser

class FragmentParser(object):

    def __init__(self, fragment_in, frag_id):

        self.atom_parsers = []

        start_index_dict = {}
        
        for atom_in in self.split_fragments_in(fragment_in):

            atom_type = "".join(itertools.takewhile(str.isupper, atom_in))
            count = int(atom_in[len(atom_type):])
            try:
                start_index = start_index_dict[atom_type]
                start_index_dict[atom_type] += count
            except KeyError:
                start_index = 1
                start_index_dict[atom_type] = 1 + count

            self.atom_parsers.append(AtomTypeParser(start_index, atom_in))

        self.frag_id = frag_id

    def get_fragment_in(self):
        return "".join([atom_parser.get_atom_in() for atom_parser in self.atom_parsers])

    def get_fragment_id(self):
        return self.frag_id

    def get_intra_molecular_variables(self):
        for atom_parser in self.get_atom_and_virtual_site_types():

            yield from [[atom1, self.get_fragment_id(), atom2, self.get_fragment_id(),
                    "x-intra-{}".format(var_type)] for atom1, atom2, var_type
                    in atom_parser.get_intra_atom_type_variables()]

        for atom_parser_index, atom_parser1 in enumerate(self.get_atom_and_virtual_site_types()):

            for atom_parser2 in list(self.get_atom_and_virtual_site_types())[atom_parser_index + 1:]:

                yield from [[atom1, self.get_fragment_id(), atom2, self.get_fragment_id(),
                        "x-intra-{}".format(var_type)] for atom1, atom2, var_type
                        in atom_parser1.get_inter_atom_type_variables(atom_parser2)]
                


    def get_inter_molecular_variables(self, other):
        for atom_parser_self in self.get_atom_and_virtual_site_types():

            for atom_parser_other in other.get_atom_and_virtual_site_types():

                yield from [[atom1, self.get_fragment_id(), atom2, other.get_fragment_id(),
                        "x-inter-{}".format(var_type)] for atom1, atom2, var_type
                        in atom_parser_self.get_inter_atom_type_variables(atom_parser_other)]

    def get_atom_and_virtual_site_types(self):

        for atom_parser in self.atom_parsers:

            yield atom_parser

    def get_atom_types(self):

        for atom_parser in self.atom_parsers:

            if atom_parser.atom_type not in constants.lone_pairs:
                yield atom_parser

    def split_fragments_in(self, fragment_in):

        # the first character (inclusive) of the current atom type
        start_index = 0

        # the last character (exclusive) of the current atom type
        end_index = 0

        # true if we are in the process of reading the atom type, false if we are in the process of reading the
        # number of atoms of a type
        reading_atom_type = True

        while(True):

            try:
                # read character code of current character
                character_code = ord(fragment_in[end_index])

            # IndexError will be thrown when we reach the end of the string
            except IndexError:

                # if we were reading atom number, then yield the last fragment and exit the generator
                if not reading_atom_type:

                    yield fragment_in[start_index:end_index]
                    return

                # if we were reading atom type, then raise an exception, because every atom type must be followed
                # by an atom number
                raise Exception

            # if character is a digit
#MRR            if character_code >= 48 and character_code < 58:
            if fragment_in[end_index].isdigit():

                # if this is the first digit of the number
                if reading_atom_type:

                    # if we didn't read at least 1 letter for atom type, raise an error
                    if(start_index == end_index):
                        raise Exception

                    # otherwise start reading in the atom number
                    reading_atom_type = False
                    
                
                end_index += 1

            # if the character is a capital letter
#MRR            elif character_code >= 65 and character_code < 91:
            elif fragment_in[end_index].isupper():

                # if this is the first digit of the atom type
                if not reading_atom_type:

                    # yield the current atom type and atom number pair
                    yield fragment_in[start_index:end_index]

                    # reset start index to start counting the next atom type and atom number pair
                    start_index = end_index

                    reading_atom_type = True

                end_index += 1

            # if character is not a digit or capital letter
            else:
                print("error", fragment_in)
                # raise exception, because fragment_in must be all numbers and capital letters.
                raise Exception

class AtomTypeParser(object):

    def __init__(self, start_index, atom_type_in):
        self.atom_type = "".join(itertools.takewhile(str.isupper, atom_type_in))
        self.start_index = start_index
        self.count = int(atom_type_in[len(self.atom_type):])

    def get_atom_in(self):
        return self.atom_type + str(self.count)

    def get_type(self):
        return self.atom_type

    def get_count(self):
        return self.count

    def get_intra_atom_type_variables(self):
        for atom_index, atom1 in enumerate(self.get_atoms()):
            for atom2 in list(self.get_atoms())[atom_index + 1:]:
                yield atom1, atom2, "{0}+{0}".format(self.atom_type)

    def get_inter_atom_type_variables(self, other):
        for atom1 in self.get_atoms():
            for atom2 in other.get_atoms():
                yield atom1, atom2, "+".join(sorted([self.atom_type,
                        other.atom_type]))

    def get_atoms(self):

        for i in range(self.start_index, self.start_index + self.count):

            yield self.atom_type + str(i)
