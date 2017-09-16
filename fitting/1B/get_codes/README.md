# Run Test
This code generates, compiles and runs a test of the 1B polynomials.
To test, just run the script run_test.sh.

The expected output of the eval code does some checks and calculates different
contributions to the energy, and numerical and analytical gradients:
 - E_nograd is the total energy obtained with the polynomials without 
   gradients
 - E_poly is the energy ONLY from the polynomials 
 - E_elec is the intramolecular electrostatic energy
 - E_disp is the intramolecular dispersion energy
 - E_grad is the total energy obtained with the polynomials with gradients
 - E_p/2p/m/2m is the energy once the x/y/z coordinate of that atom is changed
   by (+/-)1 or (+/-)2 times the step (set to 1E-6 by default)
 - Atom   Analit: XXX Numerical: YYY Diff: ZZZ 
   This line gives the analytical (from polynomials), numerical (from finite
   differences) gradient, and their absolute value of the difference.

# Run fit
 `./fit-1b [training_set.xyz] > [log_file]`
 In destination/ , this fits the polynomials to your training set in xyz format.
 - The energies of the training set must be in kcal/mol in the comment line.
 - The energies of the training set must be scaled to the minimum energy, being
   that one 0.
 - The output will be the file fit-1b.cdl. This file contains the parameters.

# Run eval
 `./eval-1b [fit-1b.nc] [config.xyz] > output`
 - You need the nc file of the training set. To obtain it:
 ```
 ncgen -o fit-1b.nc fit-1b.cdl
 ```
 - The output will contain all the energies (E_nograd, E_poly ...) previously described
