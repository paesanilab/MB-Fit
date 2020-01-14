# MB-nrg potential energy function generation

[![Build Status](https://travis-ci.org/paesanilab/potential_fitting.svg?branch=master)](https://travis-ci.org/paesanilab/potential_fitting)

This repository contains all of our codes to generate MB-nrg potential
energy functions (1B, 2B, 3B):

1. Training set generator consisting of
   - Configuration generator
   - QM MB energy calculator
2. Potential energy function generator consisting of
   - Code to generate polynomials functions
   - Fitting codes


## Changelog

Most recent changes at the top

Version v0.2.3 (Under development, not released.)
* Fixed a bug where retrieve_best_fit() would think some converged fits failed
to converge.
* Virtual site labels can now be specified by the user everwhere they are used.
Previously, some places allowed user specification, while others always used
X, Y, and Z. X, Y, and Z are still the default if the user does not specify something
else.
* MBX cpp and h files are now moved to MBX_files directory during generate_mbnrg_fitting_code()
instead of being left in the current working directory until being copied into the MBX_files
directory during fitting.generate_software_files().
* Fixed a bug where make clean would not correctly remove all executables.
* DE and alpha are now only passed in as arguments to prepare_fits(). Previously,
they were also set in the config.ini by generate_fitting_config_file_new(), which
could then be overriden by the values passed into prepare_fits(). Now, if not
passed into prepare_fits() defaults of 20 (DE) and 0.0005 (alpha) are used.
* Fixed a bug where the ratio (Effective Volume / Free Volume) would not be
raised to the power of (4/3) when calculating effective polarizabilities.

Version v0.2.2
* Fixed a bug where Qchem optimizations would crash when using the
QchemCalculator module.
* Fixed a bug where Qchem optimizations and energy calculations would disable
symmetry when using the QchemCalculator module. (This was only supposed
to happen for frequency calculations.)
* Fixed a bug where c6, d6, and A constants would be out of order in
MBX files when the mon_id argument was provided out of alphabetical order
to generate_software_files().
* Fixed a bug where molecules would be passed into polynomials in the wrong
order in MBX for 3b+.
* Fixed a bug where compilation of 3b+ fitcode would fail due to including
"dispersion.h". This include statement has been removed for 3b+.
* Fixed a bug where to_xyz() and to_ghost_xyz() tests inside the TestAtom
module would fail even though those functions were working correctly.

Version v0.2.1
<br> Initial Version, all future changes will be noted in this changelog.
