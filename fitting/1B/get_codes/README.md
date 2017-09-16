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
