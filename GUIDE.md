The Following is a Guide to work thru the entire potential fitting workflow:

*** DUE TO MY FAILURE TO LOG IN, we must use @motzu.chem.ucsd.edu instead of @kierkegaard.ucsd.edu
when sshing into the guest login for ALL STEPS BELOW. Please be more careful than I was
when typing in the username and password. 3 failed logins in a row will get EVERYONE using gecko blocked
from @motzu.chem.ucsd.edu. ***

1) Make a project directory inside projects/ containing a settings.ini file
and an unoptimized geometry. Look at the example projects for examples. If
you copy settings.ini from one of the examples, update it. You can leave most fields the same, but some must
be changed:

[files] input_geometry should be the name of your input geometry file
        poly_in_path should be a generalized geometry of your molecule. for example: 
A3.in for O3, A1B2.in of CO2, A2B4.in for N2O4, or A1B1C1.in of SCN-

[molecule]  fragments should equal the number of atoms in your monomer
            charges should equal the charge of your monomer (probably 0)
            spins should equal the spin multiplicity of your monomer (probably 1)

[config_generator] num_configs should be at least 64 (ish) (as low as 40ish may be fine,
but more configs -> better fit.) Should be on the high end of 2000 for even a mediocre fit,
but keep it smaller for tests.

Activate your psi4dev conda enviornment

1) move into the configuration_generator/norm_distribution/src and run "make"

2) Move into projects/ and call python potential_fitting.py <dir name>
where <dir name> is the name of your project directory.

3) use scp to copy the files poly-grd.maple and poly-nogrd.maple located in
<your project directory>/poly to the directory projects/<your project directory>/
in the guest login given in /network-scratch/REHS18/plab-guest-account.txt

4) ssh into the guest login give in /network-scratch/REHS18/plab-guest-account.txt

5) navigate into your project directory and run:
    maple poly-grd.maple
    maple poly-nogrd.maple

6) ssh out of the guest login with the command "exit"

7) use scp to copy the files poly-grd.c, poly-grd.cpp, poly-nogrd.c, and poly-nogrd.cpp
from the guest login to <your project directory>/poly

8) move into the poly directory and run clean-maple-c.pl located in polynomial_generation/src/ as follows:

../../../polynomial_generation/src/clean-maple-c.pl < poly-grd.c > poly-grd.cpp
../../../polynomial_generation/src/clean-maple-c.pl < poly-nogrd.c > poly-nogrd.cpp


9) copy configs.ini from projects/CO2 to projects/<your project directory> and change the following:

change [common] molecule based on your molecule.
ie A3 for O3, A1B2, for CO2, A2B4 for N2O4, etc

change [fitting] number_of_atoms appropraitely
change [fitting] number_of_electrostatic_sites = number of atoms

For the below sections indexing is based on order defined in your geometry
[fitting] excluded_pairs_12 should equal an array of all the pairs of
atoms that are bonded. ie: if atoms 0-1, and 0-2 are bonded, add [0,1] and [0,2] to the list.
[fitting] excluded_pairs_13 should equal an array of all angles in the atom
where the two atoms are at the ends of the angle. ie: if atoms 0-1-2 form an angle,
add [0,2] to the list

10) open fitting/1B/get_codes/prepare_1b_fitting_code.sh in a text editor

change lines 11 and 13 so that line 11 is the absolute path to <your project dir>/poly
and line 13 is the absolute bath to the get_codes directory, you should only have to
change the username (ebullvul to your username)

* DO NOT COMMIT THESE CHANGES *

11) use mkdir to make a fitcode directory inside <your project directory>

12) switch to fitcode conda enviornment

13) move into the fitcode directory and call

../../../fitting/1B/get_codes/prepare_1b_fitting_code.sh <abs path to .in file in your project directory> <abs path to poly directory>

14) run ./fit-1b ../training_set/training_set.xyz

15) run ncgen -o fit-1b.nc fit-1b.cdl

16) your fitcode is generated! In order to test it, talk to Marc


###################################################

Instructions to make your fitcode psi4 enviornement:

deactivate your psi4dev enviornment

update conda:
conda update conda

create enviornment:
conda create -n fit-code python=3.6 anaconda

activate:
conda activate fit-code

Your fitcode enviornment should be set up!

you can switch between fitcode enviornment and psi4dev enviornment using
conda deactivate and conda activate <enviornment name>
