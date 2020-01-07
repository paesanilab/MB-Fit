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
* Fixed a bug where some libraries needed to compile fitcode would have
the incorrect path in the Makefile on some systems.

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
