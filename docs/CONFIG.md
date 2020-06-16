# Configuration File
The configuration file is the file that contains all the necessary data to generate the fitting codes for TTM-nrg and MB-nrg. 

IMPORTANT NOTE. This file should be automatically generated from the pipeline. This document just clarifies what is each one of the entries.

## The `common` section
An example of this section is shown below:
```
[common]
molecule = A1B2_A1B2
```
The different entries are explained below:
* `molecule` is the molecule symmetry, which helps identify exchangable atoms and fragments. 

Some examples are `symmetry = ["A1B4"]` for methane monomer, `symmetry = ["A1B2","A1B2"]` for a CO2 dimer, `symmetry = ["A6B6","C1D2"]` for a benzene -- water dimer without lone pairs, and `symmetry = ["A1B2Z2","C1D2"]` for a H2O -- SO2 dimer with lone pairs. The rules are the following:

1. Symmetry names must be written in capital letters and start with A for the first atom of the first monomer. Any new atom type will be assigned the next letter of the alphabet.
2. Exchangable atoms must have the same label, even if they are in different molecules.
3. If there are virtual sites such as lone pairs that will play a role in the polynomials, they must be labels with letters X, Y, or Z.
4. If two groups inside the same molecule have the same symmetry, they should be separated. As an example, DMSO should have a symmetry `symmetry = ["A1B3_A1B3_C1D1"]`. This allows permutation within the whole methyl groups, but not within the different carbons or hydrogens individually between the two methyl groups.
5. The symmetry order MUST match the xyz order.

## The `fitting` section
An example of the fitting section is shown below:
```
number_of_atoms = 6
number_of_electrostatic_sites = 6
excluded_pairs_12 = [[[0, 1], [0, 2]], [[0, 1], [0, 2]]]
excluded_pairs_13 = [[[1, 2]], [[1, 2]]]
excluded_pairs_14 = [[], []]
charges = [[0.6837, -0.3418, -0.3418], [0.6837, -0.3418, -0.3418]]
polarizabilities = [[1.2253, 0.6808, 0.6808], [1.2253, 0.6808, 0.6808]]
polarizability_factors = [[1.2253, 0.6808, 0.6808], [1.2253, 0.6808, 0.6808]]
k_min = 0.0
k_max = 50.0
d_min = 0.0
d_max = 50.0
k_min_init = 1.0
k_max_init = 4.0
d_min_init = 1.0
d_max_init = 4.0
r_in = 7.0
r_out = 8.0
c6 = [319.9415, 221.5987, 173.3298]
d6 = [3.08949, 3.71685, 4.09252]
a = [15312.3, 20732.5, 78777.2]
var_intra = exp
var_inter = exp
var_virtual_sites = coul
alpha = 0.0005
energy_range = 20
virtual_site_labels = ['X', 'Y', 'Z']
nvars = 15
npoly = 24
```

These are the explanations of each one of the `fitting` entries:
* `number_of_atoms` is the total number of real atoms of the molecule
* `number_of_electrostatic_sites` is the total number of sites, including both real atoms and virtual electrostatic sites such as the M-site of water.
* `excluded_pairs_12` is a list with the pairs of atoms that are at one bond distance from one to another. These pairs will be excluded from non-bonded force calculations such as dispersion, buckingham and electrostatics.
* `excluded_pairs_13` is a list with the pairs of atoms that are at two bonds distance from one to another. These pairs will be excluded from non-bonded force calculations such as dispersion, buckingham and electrostatics.
* `excluded_pairs_14` is a list with the pairs of atoms that are at three bonds distance from one to another. These pairs will be excluded from non-bonded force calculations such as dispersion, buckingham and electrostatics.
* `charges` is a list of the charges of the real atoms for each monomer.
* `polarizabilities` is a list of the polarizabilities of the real atoms for each monomer.
* `polarizability_factors` is a list of the polarizability factors of the real atoms for each monomer.
* `virtual_site_labels` is a list with the characters that will represent polynomial virtual sites, such as lone pairs or midbond sites.
* `var_intra` is the variable functional form that will be used in the polynomials to describe intramonomer pairs, i.e. pairs of atoms that are in the same monomer. Options are:
    - `exp` = exp(-kr)
    - `exp0` = exp(-k(r-d0))
    - `coul` = exp(-kr)/r
    - `coul0` = exp(-k(r-d0))/r
* `var_inter` is the variable functional form that will be used in the polynomials to describe intermonomer pairs, i.e. pairs of atoms that are in different monomers. Options are:
    - `exp` = exp(-kr)
    - `exp0` = exp(-k(r-d0))
    - `coul` = exp(-kr)/r
    - `coul0` = exp(-k(r-d0))/r
* `var_virtual_sites` is the variable functional form that will be used in the polynomials to describe pairs involving lone pairs, i.e. pairs of atoms in which one of them is a virtual site used in polynomials such as lone pairs for water. Options are:
    - `exp` = exp(-kr)
    - `exp0` = exp(-k(r-d0))
    - `coul` = exp(-kr)/r
    - `coul0` = exp(-k(r-d0))/r
* `k_min` (`k_max`) is the minimum (maximum) value allowed for the constants `k` in the functional forms during the fitting procedure.
* `d_min` (`d_max`) is the minimum (maximum) value allowed for the constants `d0` in the functional forms during the fitting procedure.
* `k_min_init` (`k_max_init`) is the minimum (maximum) value allowed for the random initial guess of the `k` constant.
* `d0_min_init` (`d0_max_init`) is the minimum (maximum) value allowed for the random initial guess of the `d0` constant.
* `r_in` defines the distance at which the polynomials will start to be smoothly switched to 0.
* `r_out` defines the distance at which the polynomials are switched to 0.
* `c6` is a list of the C6 coefficients, with no repetitions, and only applicable if the number of monomers is 1 (in which case the list becomes the same as the one for the homodimer) or 2. Otherwise, it can be set to `[0.0]`. The order is expected to follow the letters in the `molecule` entry. In the case of CO2 -- CO2, as in the example above, the order is `c6 = [c6_AA, c6_AB, c6_BB]`. Note that `c6_BA` is not in there. Since the two monomers are the same, only `c6_ij` indexes `i`,`j` are taken when `j>=i`.
* `d6` is a list analogous to `c6` but containing the corresponding damping parameters for the dispersion, which is the same as the exponential factor of the Born-Mayer repulsion function `A*exp(-br)`. If not needed can be set to `[0.0]`
* `a` is a list analogous to `d6`, but containing the pre-exponential factor instead.
* `nvars` is the number of different variables in the polynomial.
* `npoly` is the number of terms in the polynomial.
* `alpha` is the ridge regression parameter that regulates how strongly the Tikhonov regularization is enforced. Larger values of alpha mean a narrower distribution of the order of magnitude of the linear constants in the polynomials.
* `energy_range` is the DE value in the expression of the weights: `w_i = [DE / (E_i - E_min + DE)]^2`. Lower values of DE will result in polynomials that are more accurate in the low energy region with the drawback of possible instabilities at higher energyes (holes). Higher values will minimize the presence of holes with the drawback of losing accuracy in the low energy region.
