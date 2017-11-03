import molecule
import dimer

'''
Reads a file as input and convert it into a list of molecules
'''

line_num = 0
dimer_input = 0
trimer_input = 0
monomer_input = 0

def readfile(datafile):
    # Open the input file manually
    inputfile = open(datafile,'r')
    count = 0
    fileline = 1
    compounds = []
    atom_count_arr = [1]

    # Determine what to look for
    if "dimer" in datafile:
        dimer_input = 1
    elif "trimer" in datafile:
        trimer_input = 1
    # By default, set monomer to true
    else:
        global monomer_input
        monomer_input = 1

    while fileline:
        # Empty molecules list
        molecules = []

        # This looks for atom amount
        fileline = read(inputfile)

        # Check for EOF; if there's no more to read, return
        if fileline == '':
            return compounds

        # If the line wasn't EOF, cast into list of int
        atom_count_arr = [int(fileline)]

        # Comment, if dimer or trimer, convert comment into list
        fileline = read(inputfile)
        if not monomer_input:
            atom_count_arr = fileline.split()
        for mol_count in range(len(atom_count_arr)):

            # Create a new molecule every time
            new_mol = form_mol(atom_count_arr[mol_count], inputfile)
        
            # Some debug into
            count += 1
            # print(new_mol)
            # print('That was molecule #{}'.format(count))

            # Append molecule into molecules array; this is size 1-3
            molecules.append(new_mol)
      
        # Whether it's a monomer, dimer, trimer, append to compounds
        if len(molecules) == 2:
            new_comp = dimer.Dimer(molecules)
            # print(new_comp)
            compounds.append(new_comp)
        else:
            compounds.append(molecules[0])
        #print(compounds)

    inputfile.close()

def read(readfile):
    """ An overriding read function that also tells line number """
    line = readfile.readline()

    global line_num
    line_num += 1
    
    return line

'''
Make molecules
'''
def form_mol(atom_count, inputfile):
    """ Form molecules based on number of atoms specified by input file
        Returns a molecule that can be by itself, or used to form a dimer
        or trimer.
    """
    
    # Prepare a list to append to molecule
    atom_arr = []

    # Read info about each atom
    for atom in range(int(atom_count)):
        fileline = read(inputfile)
        atom_info = fileline.split()

        # Form the atom and add it to list
        atom_arr.append(molecule.Atom(atom_info))  

    # Form the molecule
    new_mol = molecule.Molecule(atom_arr)

    #print(new_mol)
    return new_mol
