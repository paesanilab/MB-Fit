# Smile Strings for MBML.

## Where do I have to specify the SMILE string?

You must specify the property SMILES inside the [molecule]
section of your '.ini' settings file. The SMILE strings
for your molecule are specified as a comma delimited list
of the SMILE strings for each fragment. <br>
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

## Why do I have to specify the SMILE string?

The SMILE string for each fragment is used to
find the 12, 13, and 14 excluded pairs of the fragment. <br>
In order to find which atoms are 1, 2, or 3 bonds apart,
we must know which atoms are bonded to eachother. This 
information is included in the SMILE string. <br>
<br>
Additionally, we need to know which atoms are bonded
to apply standard order before interacting with the database,
so we need the SMILE string for this too.

## How should I specify the SMILE string?

The SMILE string for each fragment must follow the
[general rules](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system#Description)
of SMILE strings in addition to the following extra rules:
 * Order of atoms in the SMILE string must match the
 order of the atoms in the '.xyz' file.
 * All Hydrogens must be included explicity.
 * Do not include any information about charges
 in your SMILE string.
 
## Is there an algorithm I can apply to generate these SMILE strings?

Yes! While it doesn't matter how you generate the SMILE
string as long as it follows the above rules, this is the
algorithm I recommend applying:

1) Write out all atoms of the fragment in the
same order they are in your '.xyz' file. Remember, if
any atom name is 2 letters, you must but braces around it.

2) Add a '.' between any atoms that are not bonded.

3) Add numbers after atoms to signify bonds between
non-adjacent atoms. In the SMILE format you can bond two
non-adjacent atoms by putting the same digit after them.
For example, the SMILE string 'C1CC1' species a ring of 3
carbons. When there are multiple digits after a single atom
each digit is treated independently, 'C12.C23.C31' also
specifies a ring of 3 carbons. If you need to specify
more than 10 bonds this way you can use % signs to
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
Ste p2) CO.H.H
Step 3) C12O.H1.H2

Here is another example:
```
3
comment line
Ph   x   y   z
O   x   y   z
O   x   y   z
```

Step 1) [Ph]OO
Step 2) [Ph]O.O
Step 3) [Ph]1O.O1

## Any feedback?

Let Ethan know if any of the above is confusing, so I
can update this doc!
