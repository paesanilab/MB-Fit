import molecule
import dimer

'''
Reads a file as input and convert it into a list of molecules
'''

lineNum = 0
dimer = 0
trimer = 0
monomer = 0

def readfile(datafile):
    # Open the input file manually
    inputfile = open(datafile,'r')
    count = 0
    fileline = 1
    compounds = []
    atomCountArr = [1]

    # Determine what to look for
    if "dimer" in datafile:
        dimer = 1
    elif "trimer" in datafile:
        trimer = 1
    # By default, set monomer to true
    else:
        global monomer
        monomer = 1

    while fileline:
        # Empty molecules list
        molecules = []

        # This looks for atom amount
        fileline = read(inputfile)

        # Check for EOF; if there's no more to read, return
        if fileline == '':
            return molecules

        # If the line wasn't EOF, cast into list of int
        atomCountArr = [int(fileline)]

        # Comment, if dimer or trimer, convert comment into list
        fileline = read(inputfile)
        if not monomer:
            atomCountArr = fileline.split()
        for molCount in range(len(atomCountArr)):

            # Create a new molecule every time
            newMol = formMol(atomCountArr[molCount], inputfile)
        
            # Some debug into
            count += 1
            #print(newMol.toString())
            print('That was molecule #{}'.format(count))

            # Append molecule into molecules array; this is size 1-3
            molecules.append(newMol)
      
        # Whether it's a monomer, dimer, trimer, append to compounds
        if len(molecules) == 2:
            newComp = dimer(molecules)
            print(newComp.toString())
            compounds.append(newComp)
        else:
            compounds.append(molecules)
        #print(compounds)

    inputfile.close()

'''
An overriding read function that also tells line number
'''
def read(readfile):
    line = readfile.readline()

    global lineNum
    lineNum += 1
    
    return line

'''
Make molecules
'''
def formMol(atomCount, inputfile):
    
    # Prepare a list to append to molecule
    atomArr = []

    # Read info about each atom
    for atom in range(int(atomCount)):
        fileline = read(inputfile)
        atomInfo = fileline.split()

        # Form the atom and add it to list
        atomArr.append(molecule.atom(atomInfo))  

    # Form the molecule
    newMol = molecule.molecule(atomArr)

    #print(newMol.toString())
    return newMol

readfile("dimer_input.xyz")
