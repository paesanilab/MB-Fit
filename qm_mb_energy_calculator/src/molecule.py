from hashlib import sha1

class Atom(object):
    """
    A class for an atom in a fragment, stores it's charge, number of unpaired
    electrons, and coordinates
    """

    '''
    Initialize a new atom from the given information
    '''
    def __init__(self, name, x, y, z):
        # Single letter symbol corresponding to periodic table abbrehviation.
        # For example: H, O, N, Cl, He
        self.name = name
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
    Returns a string representing the information in this atom in the xyz file
    format
    '''
    def to_xyz(self):
        return "{:2} {:22.14e} {:22.14e} {:22.14e}".format(self.name, self.x, self.y, self.z)

    '''
    Returns a string representing the information in this atom in the xyz file
    format, specifying this atom as a ghost atom
    '''
    def to_ghost_xyz(self):
        return "@{:2} {:22.14e} {:22.14e} {:22.14e}".format(self.name, self.x, self.y, self.z)

class Fragment(object):
    """
    A class for a fragment of a Molecule. Contains atoms
    """
    
    '''
    Initialize a new Fragment with an empty atoms list
    '''
    def __init__(self, charge, spin_multiplicity):
        # Array of atoms in this molecule
        self.atoms = []
        # charge of this fragment
        self.charge = charge
        # spin_multiplicity multiplicity of this fragment
        if spin_multiplicity < 1:
            raise ValueError("Fragment cannot have spin_multiplicity multiplicity {}. Must be greater than or equal to 1.".format(spin_multiplicity))
        self.spin_multiplicity = spin_multiplicity

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
    Get the total charge of this fragment.
    '''
    def get_charge(self):
        return self.charge

    '''
    Get the spin_multiplicity multiplicity of this fragment
    '''
    def get_spin_multiplicity(self):
        return self.spin_multiplicity

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

    '''
    Returns a string represeting the information in this fragment in the xyz
    file format specifying the atoms in this fragment as ghost atoms
    '''
    def to_ghost_xyz(self):
        string = ""
        for atom in self.atoms:
            string += atom.to_ghost_xyz() + "\n"
        return string

    '''
    Returns a string representing the information in this fragment in the xyz
    file format in STANDARD ORDER

    STANDARD ORDER is currently defined as lexigraphical order based on the
    atoms' to_xyz() string representaiton. Should probably be changed!
    '''
    def to_standard_xyz(self):
        # Get list of atom strings
        atom_strings = []
        for atom in self.atoms:
            atom_strings.append(atom.to_xyz())
        
        # Order list of atom strings
        atom_strings.sort()

        # build string to output
        string = ""
        for atom_string in atom_strings:
            string += atom_string + "\n"
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
        self.natoms = []
        self.charges = []
        self.spin_multiplicitys = []
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
        return self.fragments[:] # this copies a list in python

    '''
    Gets a list of all the atoms in this molecule
    '''
    def get_atoms(self):
        atoms = []
        for fragment in self.fragments:
            atoms += fragment.get_atoms()
        return atoms

    '''
    Get total charge of this Molecule by summing charges of Fragments.
    '''
    def get_charge(self, fragments = None):
        if fragments == None:
            fragments = range(len(self.fragments))
        charge = 0
        for index in fragments:
            charge += self.fragments[index].get_charge()
        return charge

    '''
    Get total spin_multiplicity multiplicity of this Molecule by summing spin_multiplicity
    of Fragments.
    '''
    def get_spin_multiplicity(self, fragments = None):
        if fragments == None:
            fragments = range(len(self.fragments))
        spin_multiplicity = 1
        for index in fragments:
            spin_multiplicity += self.fragments[index].get_spin_multiplicity() - 1
        return spin_multiplicity

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
    def to_xyz(self, fragments = None, cp = False):
        # by default, use all fragments
        if fragments == None:
            fragments = range(len(self.fragments))
        string = ""
        for index in range(len(self.fragments)):
            if index in fragments:
                string += self.fragments[index].to_xyz()
            elif cp:
                string += self.fragments[index].to_ghost_xyz() 
        return string[:-1] # removes last character of string (extra newline)

    '''
    Returns a string representing the fragments of this Molecule in the xyz
    file format in STANDARD ORDER.

    STANDARD ORDER is currently defined as the lexigraphical order by the
    to_standard_xyz() string representation of each fragment
    '''
    def to_standard_xyz(self):
        # get list of fragment strings
        fragment_strings = []
        for fragment in self.fragments:
            fragment_strings.append(fragment.to_standard_xyz())
        # sort the list of fragment strings
        fragment_strings.sort()
        # build string for return
        string = ""
        for fragment_string in fragment_strings:
           string += fragment_string
        return string[:-1] # removes last character of string (extra newline)

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
    Returns a string containing the many body interaction energies, in the
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
    Generates a SHA1 hash based on our molecule.
    We are using symbols, coordinates, and charges.
    Spin multiplicity to be added later.
    '''
    def get_SHA1(self):
        
        hash_string = self.to_standard_xyz() + "\n" + str(self.get_charge()) + "\n" + str(self.get_spin_multiplicity())
        return sha1(hash_string.encode()).hexdigest()
