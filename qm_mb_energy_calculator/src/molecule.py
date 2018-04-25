class Fragment(object):
    """
    A class for a fragment of a Molecule. Accepts symbols and coordinates.
    """

    def __init__(self):
        self.natoms = 0
        self.symbols = []
        self.coordinates = []
        self.comment = None

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        s = ""
        for symbol,coord in zip(self.symbols, self.coordinates):
            s += "{:2s}".format(symbol)
            for c in coord:
                s += "{:8.4f}".format(c)
            s += "\n"
        return s
 
    def get_xyz_string(self):
        """ 
        Builds a string used for an output.
        """
        s=str(self.natoms)+"\n"
        s+=self.comment+"\n"
        s+=self.__str__()
        return s

class Molecule(object):
    """
    A class for a complete assembly of fragments.
    """
    # Goal: Convert molecule into a dictionary to convert into JSON
    def __init__(self):
        self.natoms = 0
        self.fragments = []
        self.energies = {}
        self.nmer_energies = []
        self.mb_energies = []
        #TODO: consider other attributes required by this class

    def __repr__(self):
        ret_str = ""
        for frag in self.fragments:
            ret_str += str(frag)
        return ret_str[:-1]

    def write_frag_energy(self):
        """
        Outputs energies and coordinates of fragment combinations compliant to
        the xyz input file.
        """
        ret_str = ""
        for energy in self.energies.keys():
            # Issue with using "%16.8f"%: it will fill up unused bytes with
            # blanks. Presumably null characters?
            ret_str += "E{}: {}\n".format(energy, "%.8f"%self.energies[energy])
        return ret_str

    def write_mb_energy(self, limit):
        """
        Outputs the many body interaction energies.
        """
        ret_str = ""
        for i in range(limit):
            ret_str += "V_{}B: {}\n".format(i+1,"%.8f"%self.mb_energies[i])
        return ret_str

    def mol_comb_str(self, comb):
        """ 
        Builds a string representation of a part of the molecule.
        """
        ret_str = ""
        for frag_index in comb:
            ret_str += str(self.fragments[frag_index])
        return ret_str
     
    def read_xyz(self,xyzfile):
        """
        Read molecule information from xyzfile
        xyzfile must be a file handle
        """
        try:
            self.natoms = int(xyzfile.readline())

            # Since we are still storing amount of monomer atoms in comment
            # Note: Want to write an error-catching mechanism here for
            #       when number of atoms do not match with given
            self.comment = xyzfile.readline().split()

            # Error catching mechanism
            total_from_comment = 0
            for natom in self.comment:
                total_from_comment += int(natom)
            # Should the total amount not match, raise an error
            if total_from_comment != self.natoms:
                raise ValueError
          
            for mono_natoms in self.comment:
                mono_natoms = int(mono_natoms)
                new_frag = Fragment()
                for atom in range(mono_natoms):
                    tokens = xyzfile.readline().split()
                    new_frag.symbols.append(tokens[0])
                    new_frag.coordinates.append(
                        [float(coord) for coord in tokens[1:]])
                    new_frag.natoms += 1
                self.fragments.append(new_frag)

            return 1
        except ValueError:
            return 0
