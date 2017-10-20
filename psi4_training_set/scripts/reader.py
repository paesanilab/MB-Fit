import molecule

'''
Inputs a file with the input and convert it into usable info
'''

lineNum = 0;

def readfile(file):
    # Open the input file manually
    inputfile = open(file,'r')
    count = 0
    fileline = 1
    molecules = []
    while fileline:
        
        # This looks for atom amount
        fileline = read(inputfile)

        # Check for EOF
        if fileline == '':
            return molecules

        # Create an array of atoms
        atomArr = []

        # If the line wasn't EOF, cast into int
        atomCount = int(fileline)

        # Comment, skip a line
        fileline = read(inputfile)

        # Read info about each atom
        for atom in range(atomCount):
            fileline = read(inputfile)
            atomInfo = fileline.split()

            # Form the atom and add it to list
            atomArr.append(molecule.atom(atomInfo))  

        # Form the molecule
        newMol = molecule.molecule(atomArr)

        print(newMol.toString())

        # Some debug into
        count += 1
        print('That was molecule #{}'.format(count))

        molecules.append(newMol)

    inputfile.close()

'''
An overriding read function that also tells line number
'''
def read(readfile):
    line = readfile.readline()

    global lineNum
    lineNum += 1
    
    return line

readfile('input.xyz')
