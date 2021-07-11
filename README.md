# MB-Fit: Software infrastructure for data-driven many-body potential energy functions

[![Build Status](https://travis-ci.org/paesanilab/potential_fitting.svg?branch=master)](https://travis-ci.org/paesanilab/potential_fitting)

MB-Fit is an integrated software infrastructure that enables the automated development of fully transferable, data-driven many-body potential energy functions (PEFs) for generic molecules within the TTM-nrg and MB-nrg theoretical/computational frameworks. 

MB-Fit provides a complete array of software tools to: 1) generate training and test sets for individual many-body energies, 2) set up and perform the required quantum mechanical calculations of the necessary training data, 3) optimize both linear and non-linear parameters entering the mathematical expressions for the TTM-nrg and MB-nrg PEFs, and 4) generate the associated codes that are directly exported to the MBX many-body energy/force calculator (http://paesanigroup.ucsd.edu/software/mbx.html) that enables MD simulations with the TTM-nrg and MB-nrg PEFs using LAMMPS (https://www.lammps.org) and i-PI (http://ipi-code.org).

Key references for TTM-nrg and MB-nrg PEFs:
MB-pol PEF
- V. Babin, C. Leforestier, F. Paesani, “Development of a first principles water potential with flexible monomers: Dimer potential energy surface, VRT spectrum, and second virial coefficient”, J. Chem. Theory Comput. 9, 5395 (2013). https://doi.org/10.1021/ct400863t
- V. Babin, G.R. Medders, F. Paesani, “Development of a first principles water potential with flexible monomers. II: Trimer potential energy surface, third virial coefficient, and small clusters”, J. Chem. Theory Comput. 10, 1599 (2014). https://doi.org/10.1021/ct500079y
- G.R. Medders, V. Babin, F. Paesani, “Development of a first principles water potential with flexible monomers. III: Liquid phase properties”, J. Chem. Theory Comput. 10, 2906 (2014). https://doi.org/10.1021/ct5004115
- F. Paesani, “Getting the right answers for the right reasons: Towards predictive molecular simulations of water with many-body potential energy functions”, Acc. Chem. Res. 49, 1844 (2016). https://doi.org/10.1021/acs.accounts.6b00285
- S.K. Reddy, S.C. Straight, P. Bajaj, C.H. Pham, M. Riera, D.R. Moberg, M.A. Morales, C. Knight, A.W. Götz, F. Paesani, “On the accuracy of the MB-pol many-body potential for water: Interaction energies, vibrational frequencies, and classical thermodynamic and dynamical properties from clusters to liquid water and ice”. J. Chem. Phys. 145, 194504 (2016). https://doi.org/10.1063/1.4967719

TTM-nrg and MB-nrg PEFs
- D.J. Arismendi-Arrieta, M. Riera, P. Bajaj, R. Prosmiti, F. Paesani, “The i-TTM model for ab initio-based ion-water interaction potentials. I. Halide-water potential energy functions”, J. Phys. Chem. B 120, 1822 (2016). https://doi.org/10.1021/acs.jpcb.5b09562
- M. Riera, A.W. Götz, F. Paesani, “The i-TTM model for ab initio-based ion-water interaction potentials. II. Alkali metal ion-water potential energy functions”, Phys. Chem. Chem. Phys. 18, 30334 (2016). https://doi.org/10.1039/C6CP02553F
- P. Bajaj, A.W. Götz, F. Paesani, “Toward chemical accuracy in the description of ion-water interactions through many-body representations. I. Halide-water dimer potential energy surfaces”, J. Chem. Theory Comput. 12, 2698 (2016). https://doi.org/10.1021/acs.jctc.6b00302
- M. Riera, N. Mardirossian, P. Bajaj, A.W. Götz, F. Paesani, “Toward chemical accuracy in the description of ion-water interactions through many-body representations. Alkali-water dimer potential energy surfaces”, J. Chem. Phys. 147, 161715 (2017). https://doi.org/10.1063/1.4993213
- P. Bajaj, J.O. Richardson, F. Paesani, “Ion-mediated hydrogen-bond rearrangement through tunneling in the iodide–dihydrate complex”, Nat. Chem. 11, 367 (2019). https://doi.org/10.1038/s41557-019-0220-2
- P. Bajaj, D. Zhuang, F. Paesani, “Specific ion effects on hydrogen-bond rearrangements in the halide–dihydrate complexes”, J. Phys. Chem. Lett. 10, 2823 (2019). https://doi.org/10.1021/acs.jpclett.9b00899
- D. Zhuang, M. Riera, G.K. Schenter, J.L. Fulton, F. Paesani, “Many-body effects determine the local structure of Cs+ in solution”, J. Phys. Chem. Lett. 10, 406 (2019). https://doi.org/10.1021/acs.jpclett.8b03829
- F. Paesani, P. Bajaj, M. Riera, “Chemical accuracy in modeling ion-hydration through many-body representations”, Adv. Phys. X 4, 1631212 (2019). https://doi.org/10.1080/23746149.2019.1631212
- M. Riera, E.P. Yeh, F. Paesani, “Data-driven many-body models for molecular fluids: CO2/H2O mixtures as a case study”. J. Chem. Theory Comput. 16, 2246 (2020). https://doi.org/10.1021/acs.jctc.9b01175
- M. Riera, A. Hirales, F. Paesani, “Data-driven many-body models with chemical accuracy for CH4/H2O mixtures”. J. Phys. Chem. B 124, 11207 (2020). https://doi.org/10.1021/acs.jpcb.0c08728
- E. Lambros, S. Dasgupta, E. Palos, S. Swee, J. Hu, F. Paesani, “General many-body framework for data-driven potentials with arbitrary quantum mechanical accuracy: Water as a case study”. https://doi.org/10.26434/chemrxiv.14710815.v1    
- 155.	A. Caruso, F. Paesani, “Data-driven many-body models enable a quantitative description of chloride hydration from clusters to bulk”. Under review. https://doi.org/10.26434/chemrxiv.14755449.v1



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
