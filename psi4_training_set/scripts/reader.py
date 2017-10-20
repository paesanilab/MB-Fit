import molecule

'''
Inputs a file with the input and convert it into usable info
'''

def readfile(file):
    # Open the input file manually
    inputfile = open(file,'r')
    count = 0
    fileline = 1
    while fileline:
        
        # This looks for atom amount
        fileline = inputfile.readline() # Atom count, record this

        # Check for EOF
        if fileline == '':
            break

        # Create an array of atoms
        atomArr = []

        # If the line wasn't EOF, cast into int
        atomCount = int(fileline)

        # Comment, skip a line
        fileline = inputfile.readline() # Comment, skip a line

        # Read info about each atom
        for atom in range(atomCount):
            fileline = inputfile.readline()
            atomInfo = fileline.split()

            # Form the atom and add it to list
            atomArr.append(molecule.atom(atomInfo))  

        # Form the molecule
        newMol = molecule.molecule(atomArr)

        print(newMol.toString())

        # Some debug into
        count += 1
        print('That was molecule #{}'.format(count))

    inputfile.close()

readfile('input.xyz')
