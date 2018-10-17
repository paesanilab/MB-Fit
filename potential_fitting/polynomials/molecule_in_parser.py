# external package imports
import itertools

class MoleculeInParser(object):

    def __init__(self, molecule_in):
        self.fragment_parsers = []

        frag_id = 'a'

        for fragment_in in molecule_in.split("_"):

            self.fragment_parsers.append(FragmentParser(fragment_in, frag_id))

            frag_id = chr(int(frag_id) + 1)

    def get_molecule_in(self):
        return "_".join([fragment_parser.get_fragment_in() for fragment_parser in self.fragment_parsers])

    def get_fragments(self):

        for fragment_parser in self.fragment_parsers:

            yield fragment_parser

class FragmentParser(object):

    def __init__(self, fragment_in, frag_id):

        self.atom_parsers = []

        for atom_in in self.split_fragment_in(fragment_in):

            self.atom_parsers.append(AtomTypeParser(atom_in))

        self.frag_id = frag_id

    def get_fragment_in(self):
        return "".join([atom_parser.get_atom_in() for atom_parser in self.atom_parsers])

    def get_fragment_id(self):
        return self.frag_id

    def get_intra_molecular_variables(self):
        for atom_parser in self.get_atom_types():

            yield from [[atom1, self.get_fragment_id(), atom2, self.get_fragment_id()]
                    for atom1, atom2 in atom_parser.get_intra_atom_type_variables()]

        for atom_parser_index, atom_parser1 in enumerate(self.get_atom_types()):

            for atom_parser1 in self.get_atom_types()[atom_parser_index:]:

                yield from [[atom1, self.get_fragment_id(), atom2, self.get_fragment_id()]
                        for atom1, atom2 in atom_parser1.get_intra_atom_type_variables(atom_parser2)]
                


    def get_inter_molecular_variables(self, other):
        for atom_parser_self in self.get_atom_types():

            for atom_parser_other in other.get_atom_types():

                yield from [[atom1, self.get_fragment_id(), atom2, other.get_fragment_id()]
                        for atom1, atom2 in atom_parser_self.get_intra_atom_type_variables(atom_parser_other)]
                

    def get_atom_types(self):

        for atom_parser in self.atom_parsers:

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
                character_code = int(fragment_in[start_index])

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
            if character_code >= 48 and character_code < 58:

                # if this is the first digit of the number
                if reading_atom_type:

                    # if we didn't read at least 1 letter for atom type, raise an error
                    if(start_index == end_index):
                        raise Exception

                    # otherwise start reading in the atom number
                    reading_atom_type = False
                    
                
                end_index++

            # if the character is a capital letter
            elif character_code >= 65 and character_code < 91:

                # if this is the first digit of the atom type
                if not reading_atom_type:

                    # yield the current atom type and atom number pair
                    yield fragment_in[start_index:end_index]

                    # reset start index to start counting the next atom type and atom number pair
                    start_index = end_index

                    reading_atom_type = True

                end_index++

            # if character is not a digit or capital letter
            else:

                # raise exception, because fragment_in must be all numbers and capital letters.
                raise Exception
class AtomTypeParser(object):

    def __init__(self, atom_type_in):
        count_str = "".join(itertools.takewhile(str.isupper, atom_type_in))
        self.count = int(count_str)
        self.atom_type = atom_type_in[len(count_str):]

    def get_atom_in(self):
        return self.atom_type + str(self.count)

    def get_intra_atom_type_variables(self):
        for i in range(self.count):
            for k in range(i, self.count):
                yield self.atom_type + str(i), self.atom_type + str(k)

    def get_inter_atom_type_variables(self, other):
        for i in range(self.count):
            for k in range(other.count):
                yield self.atom_type + str(i), other.atom_type + str(k)


    def get_atoms(self):

        for i in range(self.count):

            yield self.atom_type
