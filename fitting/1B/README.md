# Multiple fits (O. Hamto)

These codes fit a training set to polynomials made up of the morse-like variable exp(-k(d-d0)) with the nonlinear parameters "k" and "d" describing the intramolecular bonds of a given monomer. These parameters are varied during the fitting procedure. 

It may be necessary to fit a training set several times using different sets of initial guesses for the nonlinear parameters (k-intra, d-intra) in order to find the least possible weighted RMSD.

The following are guidelines for varying the initial guesses:

1. Choose an allowable range for the final values of the parameters "k-intra" and "d-intra"
	k < 4
	0 < d < 1.5*(equilibrium distances)

2. Perform multiple fits until ~20% of the weighted RMSD's agree and/or 500 fits are done.

3. Of the remaining 20%, choose the best performing fit in terms of the test set RMSD's.
