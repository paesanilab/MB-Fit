class Atom(object):
    """
    A class for an atom in a fragment, stores it's charge, number of unpaired
    electrons, and coordinates
    """

    '''
    Initialize a new atom from the given information
    '''
    def __init__(self, name, charge, unpaired_electrons, x, y, z):
        # Single letter symbol corresponding to periodic table abbrehviation.
        # For example: H, O, N, Cl, He
        self.name = name
        # Ionization of atom, positive values means missing electrons, negative
        # values mean extra electrons
        self.charge = charge
        # Number of unpaired electrons. An electron is unpaired if it occupies
        # an orbital without a second electron
        self.unpaired_electrons = unpaired_electrons
        # x position in angstroms
        self.x = x
        # y position in angstroms
        self.y = y
        # z position in angstoms
        self.z = z

    '''
    Get the name of this symbol as its periodic table abbreviation
    '''
    def get_name(self):
        return self.name

    '''
    Get the charge of this atom
    '''
    def get_charge(self):
        return charge

    '''
    Get the number of unpaired electrons in this atom
    '''
    def get_unpaired(self):
        return unpaired_electrons

    '''
    Returns a string representing the information in this atom in the xyz file
    format
    '''
    def to_xyz(self):
        return "{} {} {} {}".format(self.name, self.x, self.y, self.z)

class Fragment(object):
    """
    A class for a fragment of a Molecule. Contains atoms
    """
    
    '''
    Initialize a new Fragment with an empty atoms list
    '''
    def __init__(self):
        # Array of atoms in this molecule
        self.atoms = []

    """
    STILL NEEDED?
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        s = ""
        for symbol,coord in zip(self.symbols, self.coordinates):
            s += "{:2s}".format(symbol)
            for c in coord:
                s += "{:22.14e}".format(c)
            s += "\n"
        return s
    """

    '''
    Add an Atom to this Fragment.
    '''
    def add_atom(self, atom):
        self.atoms.append(atom)

    '''
    Get the total charge of this fragment by summing the charges of the atoms
    within it.
    '''
    def get_charge(self):
        charge = 0
        for atom in self.atoms:
            charge += atom.get_charge()
        return charge

    '''
    Get the total number of unpaired electrons in this fragment by summing
    the unpaired electrons of the Atoms within it
    '''
    def get_unpaired(self):
        unpaired = 0
        for unpaired in self.atoms:
            unpaired += atom.get_unpaired()
        return unpaired

    '''
    Gets the number of atoms in this fragment
    '''
    def get_num_atoms(self):
        return len(self.atoms)
 
    '''
    Returns a string representing the information in this fragment in the xyz
    file format.
    '''
    def to_xyz(self):
        """ 
        Builds a string used for an output.
        """
        string = ""
        for atom in self.atoms:
            string += atom.to_xyz() + "\n"
        return string

class Molecule(object):
    """
    A Molecule holds an array of fragments.
    """

    '''
    Construct a new Molecule.
    '''
    def __init__(self):
        # list of fragments in this molecule
        self.fragments = []
        self.energies = {}
        self.nmer_energies = []
        self.mb_energies = []
        #TODO: consider other attributes required by this class

    '''
    Add a Fragment to this Molecule.
    '''
    def add_fragment(self, fragment):
        self.fragments.append(fragment)

    '''
    Get total charge of this Molecule by summing charges of Fragments.
    '''
    def get_charge(self):
        charge = 0
        for fragment in self.fragments:
            charge += fragment.get_charge()
        return charge

    '''
    Get total unpaired electrons in this Molecule by summing unpaired electrons
    of Fragments.
    '''
    def get_unpaired(self):
        upaired = 0
        for fragment in self.fragments:
            unpaired += fragment.get_unpaired()

    '''
    Gets the number of Atoms in this Molecule
    '''
    def get_num_atoms(self):
        atoms = 0
        for fragment in self.fragments:
            atoms += fragment.get_num_atoms()
        return atoms
    

    '''
    Return a string representing the fragments of this Molecule specified by
    the indicies in the fragments parameter in the xyz file format.
    '''
    def to_xyz(self, fragments):
        string = ""
        for index in fragments:
            # newline improves readability but might break psi4 or qchem
            string += self.fragments[index].to_xyz() + "\n"
        return string[:-1]

    '''
    STILL NEED?

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
                print("Error: Fragment atoms do not add to total atoms in xyz file")
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
    '''
