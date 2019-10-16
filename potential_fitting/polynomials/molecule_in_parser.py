# external package imports
import itertools

# absolute module imports
from potential_fitting.utils import constants
from potential_fitting.exceptions import InvalidInputError

class FragmentSymmetryParser(object):

    def __init__(self, symmetry):
        self.set_sub_parsers(symmetry)

    def set_sub_parsers(self, symmetry):

        self.sub_parsers = []

        if "_" in symmetry or "(" in symmetry or ")" in symmetry:

            self.has_sub_fragments = True

            cur_fragment = ""

            num_open_parens = 0

            for char in symmetry:

                if char is '(':
                    if num_open_parens > 0:
                        cur_fragment += char
                    num_open_parens += 1

                elif char is ')':
                    num_open_parens -= 1
                    if num_open_parens > 0:
                        cur_fragment += char

                elif char is '_' and num_open_parens == 0:
                    self.sub_parsers.append(FragmentSymmetryParser(cur_fragment))
                    cur_fragment = ""

                else:
                    cur_fragment += char

            if num_open_parens != 0:
                raise InvalidInputError(
                    "{} is not a valid symmetry input. Make sure all parenthesis are closed!".format(symmetry))

            self.sub_parsers.append(FragmentSymmetryParser(cur_fragment))

        else:

            self.has_sub_fragments = False

            cur_atom = ""

            for char in symmetry:
                if not char.isdigit() and len(cur_atom) != 0 and cur_atom[-1].isdigit():
                    self.sub_parsers.append(AtomSymmetryParser(cur_atom))
                    cur_atom = char
                else:
                    cur_atom += char

            self.sub_parsers.append(AtomSymmetryParser(cur_atom))

    def get_sub_parsers(self):
        return self.sub_parsers

    def get_num_atoms(self):
        return sum([sub_parser.get_num_atoms() for sub_parser in self.get_sub_parsers()])

    def get_num_atoms_and_virtual_sites(self):
        return sum([sub_parser.get_num_atoms_and_virtual_sites() for sub_parser in self.get_sub_parsers()])

    def get_num_fragments(self):
        if self.has_sub_fragments:
            return len(self.get_sub_parsers())
        else:
            return 1

    def get_atoms(self):

        index_dict = {}

        for parser_index, parser in enumerate(self.get_sub_parsers()):
            fragment_id = FragmentSymmetryParser.get_fragment_id(parser_index)

            if not self.has_sub_fragments:
                fragment_id = ""

            for symmetry_class, atom_index, fragment_index in parser.get_atoms():

                try:
                    atom_index = index_dict[symmetry_class]
                    index_dict[symmetry_class] += 1
                except:
                    atom_index = 1
                    index_dict[symmetry_class] = 2

                yield (symmetry_class, atom_index, fragment_id + fragment_index)

    def get_variables(self):

        atoms = sorted(list(self.get_atoms()), key=lambda x: x[0])

        for index1, atom1 in enumerate(atoms):
            atom1_symmetry_class = atom1[0]
            atom1_index          = atom1[1]
            atom1_fragment_index = atom1[2]

            for index2, atom2 in enumerate(atoms[index1 + 1:]):
                atom2_symmetry_class = atom2[0]
                atom2_index          = atom2[1]
                atom2_fragment_index = atom2[2]

                similarity = 0
                for atom1_frag, atom2_frag in zip(atom1_fragment_index, atom2_fragment_index):
                    if atom1_frag != atom2_frag:
                        break
                    similarity += 1

                variable_type = "x-{}-{}+{}-{}".format("inter" if similarity == 0 else "intra",
                                                    atom1_symmetry_class, atom2_symmetry_class, similarity)

                yield (atom1_symmetry_class, atom1_index, atom1_fragment_index,
                      atom2_symmetry_class, atom2_index, atom2_fragment_index, variable_type)

    def get_intermolecular_variables(self):
        yield from ((atom1_symmetry_class, atom1_index, atom1_fragment_index, atom2_symmetry_class, atom2_index, atom2_fragment_index, variable_type)
                    for atom1_symmetry_class, atom1_index, atom1_fragment_index, atom2_symmetry_class, atom2_index, atom2_fragment_index, variable_type
                    in self.get_variables()
                    if atom1_fragment_index != atom2_fragment_index)

    def get_pairs(self, vsites=[]):
        pairs = set()
        for parser in self.get_sub_parsers():
            pairs = pairs.union(parser.get_pairs(vsites=vsites))

        pairs = pairs.union(set(self.get_intermolecular_pairs(vsites=vsites)))

        return list(pairs)

    def get_intermolecular_pairs(self, vsites=[]):
        pairs = set()
        for parser1, parser2 in itertools.combinations(self.get_sub_parsers(), 2):
            atoms1 = list(parser1.get_atoms())
            atoms2 = list(parser2.get_atoms())

            for atom1_symmetry, atom1_index, atom1_fragment_index in atoms1:
                for atom2_symmetry, atom2_index, atom2_fragment_index in atoms2:
                    if atom1_symmetry not in vsites and atom2_symmetry not in vsites:
                        pairs.add(tuple(sorted((atom1_symmetry, atom2_symmetry))))

        return list(pairs)

    def get_symmetry(self):
        if self.has_sub_fragments:

            return "_".join("(" + sub_parser.get_symmetry() + ")" if "_" in sub_parser.get_symmetry() else sub_parser.get_symmetry()
                            for sub_parser
                            in self.get_sub_parsers())

        return "".join(sub_parser.get_symmetry() for sub_parser in self.get_sub_parsers())

    @staticmethod
    def get_fragment_id(parser_index):
        return chr(ord('a') + parser_index)

class MoleculeSymmetryParser(FragmentSymmetryParser):

    def get_fragment_symmetries(self):
        if self.has_sub_fragments:
            return (sub_parser.get_symmetry() for sub_parser in self.get_sub_parsers())
        return (self.get_symmetry(),)

    def get_atoms(self):

        atoms = super(MoleculeSymmetryParser, self).get_atoms()

        for atom in atoms:
            if atom[2] is "":
                yield (atom[0], atom[1], 'a')
            else:
                yield atom

class AtomSymmetryParser(FragmentSymmetryParser):

    def __init__(self, symmetry):
        super(AtomSymmetryParser, self).__init__(symmetry)

        symmetry_class = ""
        num_atoms = ""

        for char in symmetry:
            if char.isdigit():
                num_atoms += char
            else:
                symmetry_class += char

        self.symmetry_class = symmetry_class
        self.num_atoms = int(num_atoms)

    def set_sub_parsers(self, symmetry):
        pass

    def get_sub_parsers(self):
        return None

    def get_num_atoms(self):

        if self.symmetry_class in ['X', 'Y', 'Z']:
            return 0

        return self.num_atoms

    def get_num_atoms_and_virtual_sites(self):
        return self.num_atoms

    def get_num_fragments(self):
        return 0

    def get_atoms(self):
        for atom_index in range(1, self.num_atoms + 1):
            yield (self.symmetry_class, atom_index, "")

    def get_variables(self):
        for index1 in range(self.num_atoms):
            for index2 in range(index1 + 1, self.num_atoms):
                atom1 = self.symmetry_class + str(index1)
                atom2 = self.symmetry_class + str(index2)
                atom1_frag = ""
                atom2_frag = ""
                variable_type = "x-intra-{}+{}-0".format(atom1, atom2)
                yield (atom1, atom2, atom1_frag, atom2_frag, variable_type)

    def get_pairs(self, vsites=[]):
        return []

    def get_intermolecular_pairs(self, vsites=[]):
        return []

    def get_symmetry(self):
        return self.symmetry_class + str(self.num_atoms)

class MoleculeInParser(object):

    def __init__(self, molecule_in):
        self.fragment_parsers = []
        frag_id = 'a'

        for fragment_in in self.split_molecule_in(molecule_in):
            self.fragment_parsers.append(FragmentParser(fragment_in, frag_id))

            frag_id = chr(ord(frag_id) + 1)

    def get_num_atoms(self):
        return sum(atom_type.get_num_atoms() for atom_type in self.get_fragments())

    def get_num_atoms_and_virtual_sites(self):
        return sum(atom_type.get_num_atoms_and_virtual_sites() for atom_type in self.get_fragments())

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

        self.fragment_in = fragment_in

    def get_num_atoms(self):
        return sum(atom_type.get_count() for atom_type in self.get_atom_types())

    def get_num_atoms_and_virtual_sites(self):
        return sum(atom_type.get_count() for atom_type in self.get_atom_and_virtual_site_types())

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

        fragment_in = fragment_in.replace("_", "").replace("(", "").replace(")", "")

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

    def get_original_fragment_in(self):
        return self.fragment_in

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
