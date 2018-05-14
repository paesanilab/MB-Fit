from hashlib import sha1

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
        return self.charge

    '''
    Get the number of unpaired electrons in this atom
    '''
    def get_unpaired(self):
        return self.unpaired_electrons

    '''
    Returns a string representing the information in this atom in the xyz file
    format
    '''
    def to_xyz(self):
        return "{:2} {:22.14e} {:22.14e} {:22.14e}".format(self.name, self.x, self.y, self.z)

    '''
    Gets a hash representing this atom, should be unique and reproducable.
    Does NOT DEPEND on charge, or unparied electrons

    PROBABLY UNNEEDED, we can just use SHA2

    Derek: Just use SHA1; collisions will be very, very rare. Also,
           you implemented the hash function in the wrong class.
    '''
    def get_hash(self):
        # define number hash
        # max atomic number is 119 (inclusve)
        h_number += self.get_atomic_number()

        # define x hash
        # max x is 99 (inclusive), but goes up to 199 for negative values
        h_x += 120 * (int)(self.x * 10 ** 8) # x * 10^8 (120 is hash multiplier for x)
        if self.x < 0:
            h_x = -h_x + 120 * (int)(100 * 10 ** 8) # if x was negative, change hash to positive and add offset
        
        # define y hash
        # max y is 99 (inclusive), but goes up to 199 for negative values
        h_y += 31900000000 * (int)(self.y * 10 ** 8) # x * 10^8 (120 is hash multiplier for x)
        if self.x < 0:
            h_y = -h_y + 31900000000 * (int)(100 * 10 ** 8) # if x was negative, change hash to positive and add offset
        
        # define z hash
        # max z is 99 (inclusive), but goes up to 199 for negative values
        h_z += 5180000000000000000 * (int)(self.z * 10 ** 8) # x * 10^8 (120 is hash multiplier for x)
        if self.x < 0:
            h_z = -h_z + 5180000000000000000 * (int)(100 * 10 ** 8) # if x was negative, change hash to positive and add offset

        # sum the hashes
        return h_number + h_x + h_y + h_z;

    
        

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

    '''
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
    '''

    '''
    Add an Atom to this Fragment.
    '''
    def add_atom(self, atom):
        self.atoms.append(atom)

    '''
    Gets a list of the atoms in this Fragment
    '''
    def get_atoms(self):
        return self.atoms[:] # copies a list in python


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
        # list of energies for this molecule, filled in by get_nmer_energies
        self.energies = {}
        # list of nmer_energies for this molecule, filled by get_nmer_energies
        self.nmer_energies = []
        self.mb_energies = []
        #TODO: consider other attributes required by this class

    '''
    Add a Fragment to this Molecule.
    '''
    def add_fragment(self, fragment):
        self.fragments.append(fragment)

    '''
    Gets a list of all the fragments in this molecule
    '''
    def get_fragments(self):
        return fragments[:] # this copies a list in python

    '''
    Gets a list of all the atoms in this molecule
    '''
    def get_atoms(self):
        atoms = []
        for fragment in self.fragments:
            atoms += atoms.fragment.get_fragments()
        return atoms

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
    Gets the number of Fragments in this Molecule
    '''
    def get_num_fragments(self):
        return len(self.fragments)

    '''
    Gets the number of Atoms in this Molecule
    '''
    def get_num_atoms(self):
        atoms = 0
        for fragment in self.fragments:
            atoms += fragment.get_num_atoms()
        return atoms
    

    '''
    Returns a string representing the fragments of this Molecule specified by
    the indicies in the fragments parameter in the xyz file format.
    '''
    def to_xyz(self, fragments = None):
        # by default, use all fragments
        if fragments == None:
            fragments = range(len(self.fragments))
        string = ""
        for index in fragments:
            # newline improves readability but might break psi4 or qchem
            string += self.fragments[index].to_xyz()
        return string[:-1]

    '''
    Returns a string containing indicies and energies of nbody fragment
    combinations in the format of the log file
    '''
    def log_frag_energy(self):
        string = ""
        # for each item in energies, add its combination indicies and energy
        # to the output string
        for combination in self.energies.keys():
            string += "E{}: {}\n".format(combination, "%.8f"%self.energies[combination])
        return string

    '''
    Returns a string containing the many body interaction energoes, in the
    format of the log file
    '''
    def log_mb_energy(self, limit):
        string = ""
        
        for index in range(limit):
            string += "V_{}B: {}\n".format(index + 1, "%.8f"%self.mb_energies[index])
        
        return string
    

    '''
    Clears the energies, nmer_energies, and mb_energies fields to make way for
    new calculations
    '''
    def clear(self):
        self.energies = {}
        self.nmer_energies = []
        self.mb_energies = []

    '''
    Gets a hash unique hash, any two molecules that are significantly different
    should have different hashes. Calling get_hash() multiple times on the same
    molecule should result in the same hash
    '''
    def get_hash(self):
        # define the hash
        h = 0

        # add the charge to the hash
        h += self.get_charge();

        # add the unpaired es to the hash
        h += 20*self.get_unpaired();

        # add atoms to the hash
        
        # get list of fragments
        fragments = self.get_fragments()

    def get_SHA1(self):
        """
        Generates a SHA1 hash based on our molecule.
        We are using symbols, coordinates, and charges.
        Spin multiplicity to be added later.
        """
        hash_string = self.to_xyz() + "\n" + str(self.get_charge())
        return sha1(hash_string).hexdigest()        
        

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
