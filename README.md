# MB-nrg potential energy function generation

[![Build Status](https://travis-ci.org/paesanilab/potential_fitting.svg?branch=master)](https://travis-ci.org/paesanilab/potential_fitting)

This repository contains all of our codes to generate MB-nrg potential
energy functions (1B, 2B). 

## Set up
In order to use MB-fit from anywhere, source the `sourceme.sh` file:
`source sourceme.sh`

Now you can import the `potential_fitting` library in python and access any of its fucntions.

Please look at the documentation and the examples to see how to run the whole workflow.


## Changelog

Most recent changes at the top

Version v0.2.3 (Under development, not released.)
* Refactored arguments for NormalModesConfigurationGenerator class constructor and
potential_fitting.generate_normal_mode_configuration funtion. Arguments should be easier to use now, but old calls to
the funtion will be broken.
* Cleaned up old files that were not needed
* Added functionality to add direct polynomials (with and without gradients) to MBX, avoinding the need of MAPLE
* Added updated examples for CO2 and NH4+
* Configuration file writing function has been improved
* MB-nrg and TTM-nrg fit code generators have been made independant.
* added new 'arguments' parameter to all functions that call psi4 or qchem. 'arguments'
is a dictionary of extra arguments to be passed on to the QM code.
* Fit code generation is now a method in the library, rather than an external call.
* Function generate_conig_new has been renamed to get_system_properties to be more explicative.
* Examples folder has been restructured with new examples and updated function calls. 
Example notebooks will be progressively added.
* Fixed a bug where the link generated in run_fit.sh by prepare_fits() would be the local path to the
training set file instead of a path to the symbolic link in the same directory.
* Fixed a bug where the training set generator would supply the wrong number of arguments
to the exponential functions when using d0 nl parameters.
* Added the Grapher class, with some functions for graphing energy distributions of TrainingSets
as well as making correlation / error graphs for TrainingSets.
* Added the Evaluator class, whose eval() method can be used to get a TrainingSet from an eval-xb
script and a training set file.
* Added the TrainingSet class, which contains a list of TrainingSetElements and functions
for performing operations on them.
* Added the TrainingSetElement class, which contains a Molecule and a dictionary of energies
and functions for performing operations on itself.
* Updated read_jobs() and JobHandler.read_all_jobs() to also have an overwrite argument
which is passed to Database.set_properties().
* Added the overwrite keyword argument to Database.set_properties() and changed it to
have the following behavior: If overwrite is False, then trying to set the energy
of a calculation which already has an energy will do nothing but give no error. If
overwrite is True, doing so will replace the existing energy. Previous behavior was
to give an error message.
* Exceptions raised from the Database class due to postgres errors no longer print
the context, as this is of no use to a user who does not know how the database works.
* Fixed a bug that would cause generate_config_file_new() to loop infinitely when
SCF fails to converge.
* Fixed a bug where delete_calculations() and delete_all_calculations() would not
work when delete_complete_calculations was set to True.
* correlation.dat files produced by execute_fits() now have a header to tell
you what is in each column.
* Fixed a bug where retrieve_best_fit() would think some converged fits failed
to converge.
* Charges, Polarizabilities, and C6 constants generated in the config.ini file
by generate_fitting_config_file_new() now round to 4 digits after the decimal
by default. generate_fitting_config_file_new() has a new argument called num_digits
that can be set to override the default value of 4.
* Added to_string() and to_string_only_value() functions to the DistributionFunction
class and all subclasses.
* All ConfigurationGenerators now print information about the temp, A, and/or
distance distributions they use.
* Fixed a bug where some libraries needed to compile fitcode would have
the incorrect path in the Makefile on some systems.
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
