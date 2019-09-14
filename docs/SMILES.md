# SMILES Strings for MBML.

## Where do I have to specify the SMILES string?

You must specify the property SMILES inside the [molecule]
section of your '.ini' settings file. The SMILES string
for your molecule is specified as a comma delimited list
of the SMILES string for each fragment. <br>
For example, for a water monomer:
```
[molecule]
SMILES = O(H)H
```
For a water dimer:
```
[molecule]
SMILES = O(H)H,O(H)H
```

## Why do I have to specify the SMILES string?

The SMILES string for each fragment is used to
find the 12, 13, and 14 excluded pairs of the fragment. <br>
In order to determine which atoms are separated by 1, 2, or 3 bonds,
we must know which atoms are bonded to each other. This 
information is contained in the SMILES string. <br>
<br>
Additionally, we need to know which atoms are bonded
to apply standard order before interacting with the database
to ensure that atom coordinates are retrieved in the correct order.
This standard order is also deduced from the molecular connectivity
that is encoded in the SMILES string.

## How should I specify the SMILES string?

The SMILES string for each fragment must follow the
[general rules](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Description)
of SMILES strings in addition to the following extra rules:
 * Order of atoms in the SMILES string must match the
 order of the atoms in the coordinate specification ('.xyz' file).
 * All hydrogen atoms must be included explicity.
 * Do not include any information about charges
 in your SMILES string.
 
## Is there an algorithm I can apply to generate these SMILES strings?

Yes! While it does not matter how you generate the SMILES
string as long as it follows the above rules, this is the
recommended algorithm:

1) Write out all atoms of the fragment in the
same order they appear in your '.xyz' file. Remember, if
any atom name is 2 letters, you must put square brackets around it.

2) Add a '.' between any atoms that are not bonded.

3) Add numbers after atoms to signify bonds between
non-adjacent atoms. In the SMILES format you can bond two
non-adjacent atoms by putting the same digit after them.
For example, the SMILES string 'C1CC1' species a ring of 3
carbons. When there are multiple digits after a single atom
each digit is treated independently. For instance, 'C12.C23.C31'
also specifies a ring of 3 carbons. In this representation we 
have removed the implicit bond between adjacent atoms by adding 
a dot, and then introduced bonds between each pair of carbon atoms.
If you need to specify more than 10 bonds you can use % signs to
combine multiple digits: 'C%10%20.C%20%30.C%30%10' also
specifies a ring of 3 carbons.

Here is an example:

Lets say our '.xyz' file for formaldehyde looks like this:
```
4
comment line
C   x   y   z
O   x   y   z
H   x   y   z
H   x   y   z
```

Step 1) COHH
Step 2) CO.H.H
Step 3) C12O.H1.H2

Here is another example:
```
3
Chloroform
C   x   y   z
H   x   y   z
Cl  x   y   z
Cl  x   y   z
Cl  x   y   z
```

Step 1) CH[Cl][Cl][Cl]
Step 2) CH.[Cl].[Cl].[Cl]
Step 3) C123H.[Cl]1.[Cl]2.[Cl]3

## Any feedback?

Let Ethan know if any of the above is confusing, so we
can update this doc!
