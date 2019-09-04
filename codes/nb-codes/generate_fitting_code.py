import sys, os, json
from potential_fitting.utils import SettingsReader
from potential_fitting.utils import constants
import utils_nb_fitting
import file_writter_nb_fitting


if len(sys.argv) != 6:
    print("Usage: ./script <settings.ini> <config.ini> <poly-direct.cpp_with_path> <degree> <poly.in name>")
    sys.exit()
else:
    settings_path = sys.argv[1]
    config_path = sys.argv[2]
    directcpp = sys.argv[3]
    degree = int(sys.argv[4])
    name = sys.argv[5]

# In[ ]:

settings = SettingsReader(settings_path)
config = SettingsReader(config_path)

################################################################################
## MONOMER PROPERTIES ##########################################################
################################################################################

# Obtain monomers from settings
monomers = settings.get("molecule", "symmetry").split(",")
number_of_monomers = len(monomers)

# Store number of atoms in each monomer in nat1 and nat2
number_of_atoms = [int(i) for i in settings.get("molecule", "fragments").split(",")]

# Set the number of sites
number_of_sites = number_of_atoms

# Get which monomer should be mb-pol monomer (if any)
use_mbpol = [int(i) for i in settings.get("molecule", "use_mbpol").split(",")]

# Define if lone pairs are used based on monomer names
# Update number of sites if needed
use_lonepairs = [0]*number_of_monomers
for i in range(number_of_monomers):
    if "X" in monomers[i] or "Y" in monomers[i] or "Z" in monomers[i]:
        use_lonepairs[i] = 1
    if use_mbpol[i] != 0:
        number_of_sites[i] += 1
    
# Obtain the lists with the excluded pairs
excluded_pairs_12 = config.getlist("fitting", "excluded_pairs_12")
excluded_pairs_13 = config.getlist("fitting", "excluded_pairs_13")
excluded_pairs_14 = config.getlist("fitting", "excluded_pairs_14")

#Obtain charges (in the order of input), pols and polfacs
charges = config.getlist("fitting", "charges")
polarizabilities = config.getlist("fitting", "polarizabilities")
polarizability_factors = config.getlist("fitting", "polarizability_factors")

################################################################################
## DIMER   PROPERTIES ##########################################################
################################################################################

# Obtain A, C6 and d6 from user in the same order as the given pairs AA, AB ...
# By default, if not specified in the config.ini file, these lists
# will be set to a list of the same length, but with all coefficients 0.0
# __FOR__ETHAN__ 
#c6_constants = config.getlist("fitting", "c6")
#C6 = c6_constants[len(c6_constants) - 1]
C6 = config.getlist("fitting", "c6")

# last element of d6_constants is list of the inter_molecular d6 constants
#d6_constants = config.getlist("fitting", "d6")
#d6 = d6_constants[len(d6_constants) - 1]

# last element of A_constants is list of the inter_molecular A constants
# Initialize A and b to 0, for now
try:
    A_buck = config.getlist("fitting", "A")
    b_buck = config.getlist("fitting", "d6")
    d6 = b_buck
except:
    A_buck = [0.0] * len(C6)
    b_buck = [0.0] * len(C6)
    d6 = [0.0] * len(C6)
#A_constants = config.getlist("fitting", "A")
#Abuck = A_constants[len(A_constants) - 1]


################################################################################
## POLYNOMIAL PROPERTIES #######################################################
################################################################################

#Ask the user for the max value of k and d
# These are the non-linear parameters of the polynomials
k_min = config.get("fitting", "k_min")
k_max = config.get("fitting", "k_max")
k_min_intra = k_min
k_max_intra = k_max

d_min = config.get("fitting", "d_min")
d_max = config.get("fitting", "d_max")
d_min_intra = d_min
d_max_intra = d_max

# Obtain inner and outer cutoff from config.ini
# This must be required 
ri = 7.0 # polynomials start to decrease
ro = 8.0 # polynomials are completely removed

# Define list of variables that are fictitious
virtual_sites_poly = config.getlist("fitting", "virtual_site_labels")

# Define kind of variables for intra, inter and lone pairs
# Options are:
# exp e^-kd
# exp0 e^-k(d-d0)
# coul0 [e^-k(d-d0)]/r
# coul [e^-kd]/r
# Recomendation is to use exp for intra and inter and coul for lone pairs
var_intra = config.get("fitting", "var_intra")
var_virtual_sites = config.get("fitting", "var_virtual_sites")
var_inter = config.get("fitting", "var_inter")

npoly = config.getint("fitting", "npoly")

# Define Energy Range for the fitting
E_range = config.getfloat("fitting", "energy_range")


################################################################################
## Prepare pair information ####################################################
################################################################################

# Obtain types for each monomer
monomer_atom_types = []
for mon in monomers:
    monomer_atom_types.append(utils_nb_fitting.get_atom_types(mon)) 

# create dictionary mapping from atom index to atom name
atom_list = []

for i in range(len(monomers)):
    atom_list.append([])
    for type_index in range(0, len(monomer_atom_types[i]), 2):
        for atom_index in range(1, int(monomer_atom_types[i][type_index + 1]) + 1):
            atom_list[i].append(monomer_atom_types[i][type_index] + str(atom_index))
    
# read the poly.in file to get the variables
variables, intra_poly_pairs, inter_poly_pairs = utils_nb_fitting.read_poly_in(name, virtual_sites_poly, var_intra, var_inter, var_virtual_sites)
nvars = len(variables)

################################################################################
## Prepare non-linear terms information ########################################
################################################################################

# Get the different non-linear parameters
nlparam_intra, nlparam_inter, nl_param_all = utils_nb_fitting.get_non_linear_parameters(variables)
# Save number in num_nonlinear
num_nonlinear = len(nl_param_all)


################################################################################
## Write monomer classes header and cpp ########################################
################################################################################

# Write the header files
for i in range(number_of_monomers):
    file_writter_nb_fitting.write_monomer_class_header(i+1)

# Write the cpp files
for i in range(number_of_monomers):
    file_writter_nb_fitting.write_monomer_class_cpp(i+1,number_of_sites[i],number_of_atoms[i],excluded_pairs_12[i],excluded_pairs_13[i],excluded_pairs_14[i], charges[i], polarizabilities[i],polarizability_factors[i])

# Write mb-pol water monomer if applicable
for i in range(number_of_monomers):
    if use_mbpol[i] == 1:
        mon_index = i+1
        file_writter_nb_fitting.write_mbpol_monomer(mon_index)


################################################################################
## Write polynomial holders header and cpp #####################################
################################################################################

system_name = monomers[0]
for i in range(1,number_of_monomers):
    system_name += "_" + monomers[i]

file_writter_nb_fitting.write_fit_polynomial_holder_header(system_name, number_of_monomers, nl_param_all, ri, ro)

file_writter_nb_fitting.write_fit_polynomial_holder_cpp(system_name, monomer_atom_types, number_of_monomers, number_of_atoms, virtual_sites_poly, use_lonepairs, nl_param_all, variables, nvars, ri, ro, k_min_intra, k_max_intra, k_min, k_max, d_min_intra, d_max_intra, d_min, d_max)


################################################################################
## Write dispersion header and cpp #############################################
################################################################################

file_writter_nb_fitting.write_dispersion_header(monomer_atom_types, virtual_sites_poly, C6, d6)

file_writter_nb_fitting.write_dispersion_cpp(monomer_atom_types, virtual_sites_poly, excluded_pairs_12[0], excluded_pairs_13[0],excluded_pairs_14[0])


################################################################################
## Write buckingham header and cpp #############################################
################################################################################
file_writter_nb_fitting.write_buckingham_header(monomer_atom_types, virtual_sites_poly, A_buck, b_buck)

ff = open("buckingham.cpp",'w')
ff.write("hello\n")
ff.close()
file_writter_nb_fitting.write_buckingham_cpp(monomer_atom_types, virtual_sites_poly, excluded_pairs_12[0], excluded_pairs_13[0],excluded_pairs_14[0])

################################################################################
## Polynomial header and cpp ###################################################
################################################################################

file_writter_nb_fitting.write_poly_fit_header(number_of_monomers, system_name, degree, nvars, npoly)

file_writter_nb_fitting.write_poly_fit_cpp(number_of_monomers, system_name, nvars, npoly, directcpp)

################################################################################
## Fitting code ################################################################
################################################################################

file_writter_nb_fitting.write_fitting_code(number_of_monomers, number_of_atoms, number_of_sites, system_name, nl_param_all, k_min, k_max, d_min, d_max, k_min_intra, k_max_intra, d_min_intra, d_max_intra)

################################################################################
## Evaluation code #############################################################
################################################################################

file_writter_nb_fitting.write_eval_code(number_of_monomers, number_of_atoms, number_of_sites, system_name)

################################################################################
## TTM-nrg fitting code ########################################################
################################################################################
if number_of_monomers == 2:
    file_writter_nb_fitting.write_fitting_ttm_code(monomer_atom_types, virtual_sites_poly, number_of_monomers, number_of_atoms, number_of_sites, system_name, k_min, k_max)

################################################################################
## TTM-nrg evaluation code #####################################################
################################################################################

    file_writter_nb_fitting.write_eval_ttm_code(monomer_atom_types, virtual_sites_poly, number_of_monomers, number_of_atoms, number_of_sites, system_name)

################################################################################
## Makefile ####################################################################
################################################################################

file_writter_nb_fitting.write_makefile(number_of_monomers, system_name)


#
#
## ## Buckingham.cpp
#
## In[ ]:
#
#
#cppname = "buckingham.cpp"
#ff = open(cppname,'w')
#a = """
##include "buckingham.h"
#
#x2b_buck::x2b_buck() {
#  xyz1 = new double[3];
#  xyz2 = new double[3];
#}
#x2b_buck::~x2b_buck() {
#  delete[] xyz1;
#  delete[] xyz2;
#}
#
#x2b_buck::x2b_buck(double * c1, double * c2, size_t n1, size_t n2) {
#  xyz1 = new double[3*n1];
#  xyz2 = new double[3*n2];
#  std::copy(c1, c1 + 3*n1, xyz1);
#  std::copy(c2, c2 + 3*n2, xyz2);
#}
#
#double x2b_buck::get_buckingham() {
#
#  double ebuck = 0.0;
#"""
#ff.write(a)
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_a),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_a[i+1])):
#        if not types_a[i] in vsites:
#            ff.write('    const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xyz1 + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_b),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_b[i+1])):
#        if not types_b[i] in vsites:
#            ff.write('    const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xyz2 + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_a),2):
#    na = 1
#    for j in range(int(types_a[i+1])):
#        for k in range(0,len(types_b),2):
#            nb = 1
#            for l in range(int(types_b[k+1])):
#                
#                if types_a[i] not in vsites and types_b[k] not in vsites:
#                    t = "".join(sorted([types_a[i], types_b[k]]))
#
#                    ff.write('  ebuck += buck(m_A_' + t  + ', m_b_' + t + ', ' + types_a[i] + "_" + str(na) + "_a" + ', ' + types_b[k] + "_" + str(nb) + "_b" + ');\n')
#
#                nb += 1
#        na += 1
#
#a = """
#  return ebuck;
#}
#
#double x2b_buck::get_buckingham(double * grd) {
#
#  double ebuck = 0.0;
#"""
#ff.write(a)
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_a),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_a[i+1])):
#        if not types_a[i] in vsites:
#            ff.write('    const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xyz1 + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_b),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_b[i+1])):
#        if not types_b[i] in vsites:
#            ff.write('    const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xyz2 + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_a),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_a[i+1])):
#        if not types_a[i] in vsites:
#            ff.write('    double* ' + types_a[i] + '_' + str(n) + '_a_g' + '= grd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
## loops over each type of atom in the input
#for i in range(0,len(types_b),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_b[i+1])):
#        if not types_b[i] in vsites:
#            ff.write('    double* ' + types_b[i] + '_' + str(n) + '_b_g' + '= grd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_a),2):
#    na = 1
#    for j in range(int(types_a[i+1])):
#        for k in range(0,len(types_b),2):
#            nb = 1
#            for l in range(int(types_b[k+1])):
#                
#                if types_a[i] not in vsites and types_b[k] not in vsites:
#
#                    t = "".join(sorted([types_a[i], types_b[k]]))
#
#                    ff.write('  ebuck += buck(m_A_' + t  + ', m_b_' + t + ', ' + types_a[i] + "_" + str(na) + "_a" + ', ' + types_b[k] + "_" + str(nb) + "_b," + types_a[i] + "_" + str(na) + "_a_g," + types_b[k] + "_" + str(nb) + "_b_g" + ');\n')
#
#                nb += 1
#        na += 1
#
#a = """
#  return ebuck;
#}
#"""
#ff.write(a)
#ff.close()
#
#
#
#
## ## Fitting routine




## ## Fitting with TTM as underlying
#
## In[ ]:
#
#
#cppname = "fit-2b-wbuck.cpp"
#ff = open(cppname,'w')
#a = """
##include <cmath>
##include <cassert>
##include <cstdlib>
##include <ctime>
#
##include <fstream>
##include <sstream>
##include <iomanip>
##include <iostream>
##include <stdexcept>
##include <chrono>
#
##include <gsl/gsl_multimin.h>
#
##define RIDGE_REGRESSION 1
#
##ifdef RIDGE_REGRESSION
##   include "rwlsq.h"
##else
##   include "wlsq.h"
##endif
#
##include "mon1.h"
##include "mon2.h"
##include "fit-utils.h"
##include "training_set.h"
##include "x2b_""" + mon1 + """_""" + mon2 + """_v1.h"
##include "electrostatics.h"
##include "coulomb.h"
##include "dispersion.h"
##include "buckingham.h"
#
#
##define dont_be_VERBOSE yes
#
#//#define DEBUG 
#
##ifdef DEBUG
##define PR(x) std::cout << #x << ": " << (x) << std::endl;
##else
##define PR(x)
##endif /* DEBUG */
#
#using namespace std;
#
#namespace {
#
#static x2b_""" + mon1 + """_""" + mon2 + """::x2b_""" + mon1 + """_""" + mon2 + """_v1x model;
#
#
##ifdef RIDGE_REGRESSION
#const double alpha = 0.0005;
##endif
#
#// ##DEFINE HERE## energy range
#const double E_range = """ + str(E_range) + """; // kcal/mol
#
#//----------------------------------------------------------------------------//
#
#static std::vector<tset::dimer> training_set;
#static std::vector<double> elec_e;
#static std::vector<double> disp_e;
#static std::vector<double> buck_e;
#static double* ts_weights = 0;
#
#namespace linear {
#
#//----------------------------------------------------------------------------//
#
#static double* A(0);
#static double* y(0);
#static double* params(0);
#
#//----------------------------------------------------------------------------//
#
#void allocate()
#{
#    A = new double[training_set.size()*model.nparams()];
#    y = new double[training_set.size()];
#
#    params = new double[model.nparams()];
#}
#
#//----------------------------------------------------------------------------//
#
#double compute_chisq(const gsl_vector* X, void* unused)
#{
#    model.set_nonlinear_parameters(X->data);
#
#    if (model.nonlinear_parameters_out_of_range())
#        return 1.0e+6;
#
#    PR(model_type::nparams())
#    for (size_t n = 0; n < training_set.size(); ++n) {
#        double mmm[model.nparams()];
#        y[n] = training_set[n].energy_twobody
#               - model.basis(training_set[n].xyz, mmm);
#        for (size_t p = 0; p < model.nparams(); ++p) {
#            A[p + n*model.nparams()] = mmm[p];
#        }
#    }
#
##   ifdef VERBOSE
#    std::cout << "=== calling wlsq::solve() ["
#              << kit::wlsq::implementation()
#              << "] ===" << std::endl;
##   endif
#
#    double chisq;
#
#//    try {
##     ifdef RIDGE_REGRESSION
#      double penaltysq;
#      kit::rwlsq::solve(training_set.size(), model.nparams(),
#                        A, y, ts_weights, alpha, params, chisq, penaltysq);
#
#      std::cout << "<#> chisq = " << chisq
#                << " : penaltysq = " << penaltysq
#                << std::endl;
##     else
#      int rank;
#      kit::wlsq::solve(training_set.size(), model_type::nparams(),
#                       A, y, ts_weights, params, chisq, rank);
#      std::cout << "<#> chisq = " << chisq
#                << " : rank = " << rank
#                << std::endl;
##     endif
#
#      if (!gsl_finite (chisq)) {
#        return 10000000.0;
#      }
#//    } catch (const std::exception& e) {
#//      return 1000000.0;
#//    }
##   ifdef VERBOSE
#    std::cout << "\\n--> chisq = " << chisq
#              << "\\n-->  rank = " << rank
#              << '\\n' << std::endl;
##   endif
#
#    return chisq;
#}
#
#//----------------------------------------------------------------------------//
#
#} // namespace linear
#
#////////////////////////////////////////////////////////////////////////////////
#
#} // namespace
#
#////////////////////////////////////////////////////////////////////////////////
#
#int main(int argc, char** argv) {
#    if (argc < 2) {
#        std::cerr << "usage: fit-2b ts1 ..."
#                  << std::endl;
#        return 0;
#    }
#      
#    long long int duration = std::chrono::duration_cast<std::chrono::milliseconds>(
#                std::chrono::system_clock::now().time_since_epoch()).count();
#
#    srand(duration);
#
#    double x0[""" + str(len(nlparam_all)) + """];
#"""
#ff.write(a)
#for i in range(len(nlparam_all)):
#    if nlparam_all[i].startswith('d'):
#        if nlparam_all[i].startswith('d_intra'):
#            ff.write('      x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(d_max_intra) - float(d_min_intra)) + ' + ' + d_min_intra + ';\n')
#        else:
#            ff.write('      x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(d_max) - float(d_min)) + ' + ' + d_min + ';\n')
#    else:
#        if nlparam_all[i].startswith('k_intra'):
#            ff.write('      x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(k_max_intra) - float(k_min_intra)) + ' + ' + k_min_intra + ';\n')
#        else:
#            ff.write('      x0[' + str(i) + '] = ((double) rand() / (RAND_MAX)) * ' + str(float(k_max) - float(k_min)) + ' + ' + k_min + ';\n')
#
#            
#a = """
#      
#
#    model.set_nonlinear_parameters(x0);
#
#    #   ifdef RIDGE_REGRESSION
#    std::cout << "<> using ridge regression with alpha = "
#              << alpha << std::endl;
##   endif
#
#    std::cout << "\\n<><><> model type = '" << model.name() << "'\\n";
#
#    {
#        const char fn[] = "fit-2b-initial.cdl";
#        std::ofstream ofs(fn);
#        std::cout << "\\n>> dumping initial model as '" << fn << "' >>\\n\\n";
#        model.as_cdl(ofs);
#    }
#
#    ++argv;
#
#    try {
#        while (--argc != 0) {
#            size_t nd = tset::load_dimers(*argv, training_set);
#            std::cout << "'" << *(argv++) << "' : "
#                      << nd << " dimers" << std::endl;
#        }
#    } catch (const std::exception& e) {
#        std::cerr << " ** Error ** : " << e.what() << std::endl;
#        return 1;
#    }
#
#    //
#    // assign weights
#    //
#
#    ts_weights = new double[training_set.size()];
#
#    double E_min;
#    size_t N_eff;
#
#    tset::setup_weights(training_set, E_range, E_min,
#                       ts_weights, N_eff, true);
#
#    std::cout << "\\n>>   E_min = " << E_min << " kcal/mol"
#                 "\\n>> E_range = " << E_range << " kcal/mol"
#                 "\\n\\n>> training set size = " << training_set.size()
#              << "\\n>>    effective size = " << N_eff << '\\n'
#              << std::endl;
#
#
#    // electrostatics
#    double tb_ref[training_set.size()];
#    for (size_t n = 0; n < training_set.size(); n++) {
#      // Saving reference 2b energies  
#      tb_ref[n] = training_set[n].energy_twobody ;
#
#      x::mon1 m1(training_set[n].xyz);
#      x::mon2 m2(training_set[n].xyz + 3*m1.get_realsites());
#
#      int system_nsites = m1.get_nsites() + m2.get_nsites();
#      int * system_is_w;
#      double* system_sitecrds;
#      double* system_charge;
#      double* system_polfac;
#      double* system_pol;
#
#      system_sitecrds = new double [system_nsites*3];
#      system_charge = new double [system_nsites];
#      system_polfac = new double [system_nsites];
#      system_pol = new double [system_nsites]; //allocates memory to the pointers
#      system_is_w = new int[system_nsites];
#      
#      std::fill(system_is_w, system_is_w + m1.get_nsites(), m1.is_w);
#      std::fill(system_is_w + m1.get_nsites(), system_is_w + system_nsites, m2.is_w);
#      
#      int * is_w_a = system_is_w;
#      int * is_w_b = system_is_w + m1.get_nsites();
#
#      std::copy(m1.get_sitecrds(), m1.get_sitecrds() + 3 * m1.get_nsites(),
#                system_sitecrds);
#      std::copy(m2.get_sitecrds(), m2.get_sitecrds() + 3 * m2.get_nsites(),
#                system_sitecrds + 3 * m1.get_nsites());
#
#      std::copy(m1.get_charges(), m1.get_charges() + m1.get_nsites(),
#                system_charge);
#      std::copy(m2.get_charges(), m2.get_charges() + m2.get_nsites(),
#                system_charge + m1.get_nsites());
#
#      std::copy(m1.get_pol(), m1.get_pol() + m1.get_nsites(),
#                system_pol);
#      std::copy(m2.get_pol(), m2.get_pol() + m2.get_nsites(),
#                system_pol + m1.get_nsites());
#
#      std::copy(m1.get_polfacs(), m1.get_polfacs() + m1.get_nsites(),
#                system_polfac);
#      std::copy(m2.get_polfacs(), m2.get_polfacs() + m2.get_nsites(),
#                system_polfac + m1.get_nsites());
#
#
#      excluded_set_type exclude12;
#      excluded_set_type exclude13;
#      excluded_set_type exclude14;
#      
#      excluded_set_type exclude12_a;
#      excluded_set_type exclude13_a;
#      excluded_set_type exclude14_a;
#      
#      excluded_set_type exclude12_b;
#      excluded_set_type exclude13_b;
#      excluded_set_type exclude14_b;
#
#      for (auto i = m1.get_begin_12(); i != m1.get_end_12(); i++) {
#        exclude12.insert(*i);
#        exclude12_a.insert(*i);
#      }
#      for (auto i = m2.get_begin_12(); i != m2.get_end_12(); i++) {
#        std::pair<size_t,size_t> p =
#                      std::make_pair(i->first + m1.get_nsites() ,
#                                     i->second + m1.get_nsites());
#        exclude12.insert(p);
#        exclude12_b.insert(*i);
#      }
#
#      for (auto i = m1.get_begin_13(); i != m1.get_end_13(); i++) {
#        exclude13.insert(*i);
#        exclude13_a.insert(*i);
#      }
#      for (auto i = m2.get_begin_13(); i != m2.get_end_13(); i++) {
#        std::pair<size_t,size_t> p =
#                      std::make_pair(i->first + m1.get_nsites() ,
#                                     i->second + m1.get_nsites());
#        exclude13.insert(p);
#        exclude13_b.insert(*i);
#      }
#
#      for (auto i = m1.get_begin_14(); i != m1.get_end_14(); i++) {
#        exclude14.insert(*i);
#        exclude14_a.insert(*i);
#      }
#      for (auto i = m2.get_begin_14(); i != m2.get_end_14(); i++) {
#        std::pair<size_t,size_t> p =
#                      std::make_pair(i->first + m1.get_nsites() ,
#                                     i->second + m1.get_nsites());
#        exclude14.insert(p);
#        exclude14_b.insert(*i);
#      }
#
#      ttm::electrostatics m_electrostatics;
#
#      ttm::smear_ttm4x smr; 
#      smr.m_aDD_intra_12 = 0.3;
#      smr.m_aDD_intra_13 = 0.3;
#      smr.m_aDD_intra_14 = 0.055;
#
#      double ener = m_electrostatics(system_nsites, system_charge, system_polfac, system_pol,
#                                    system_sitecrds, exclude12, exclude13, exclude14, 
#                                    system_is_w, smr, 0);
#      double ener_a = m_electrostatics(m1.get_nsites(), m1.get_charges(), m1.get_polfacs(), m1.get_pol(),
#                                    m1.get_sitecrds(), exclude12_a, exclude13_a, exclude14_a, 
#                                    is_w_a, smr, 0);
#      double ener_b = m_electrostatics(m2.get_nsites(), m2.get_charges(), m2.get_polfacs(), m2.get_pol(),
#                                    m2.get_sitecrds(), exclude12_b, exclude13_b, exclude14_b, 
#                                    is_w_b, smr, 0);
#
#      // Take out electrostatic energy:
#      elec_e.push_back(ener - ener_a - ener_b);
#      training_set[n].energy_twobody -= elec_e[n];
#      // std::cerr << "Conf " << n << " : Elec= " << ener ;
#
#      // Now need to take out dispersion
#      x2b_disp disp(m1.get_sitecrds(), m2.get_sitecrds(), m1.get_realsites(), m2.get_realsites());
#      ener = disp.get_dispersion();
#      disp_e.push_back(ener);
#      training_set[n].energy_twobody -= ener ;
#      
#      // And the buckingham
#      x2b_buck buck(m1.get_sitecrds(), m2.get_sitecrds(), m1.get_realsites(), m2.get_realsites());
#      ener = buck.get_buckingham();
#      buck_e.push_back(ener);
#      training_set[n].energy_twobody -= ener ;
#
#      delete[] system_sitecrds ;
#      delete[] system_charge  ;
#      delete[] system_polfac  ;
#      delete[] system_pol ;
#      delete[] system_is_w;
#
#    }
#
#        linear::allocate();
#
#    gsl_vector* x = gsl_vector_alloc(model.get_num_nonlinear_params());
#    model.get_nonlinear_parameters(x->data);
#
#    gsl_multimin_function chisq_func;
#
#    chisq_func.n = model.get_num_nonlinear_params();
#    chisq_func.f = linear::compute_chisq;
#
#    gsl_vector* ss = gsl_vector_alloc(model.get_num_nonlinear_params());
#    gsl_vector_set_all(ss, 0.1);
#
#    std::cout << "\\n<> initial simplex sides:\\n";
#    for (size_t n = 0; n < model.get_num_nonlinear_params(); ++n)
#        std::cout << n << " : " << ss->data[n] << "\\n";
#
#    gsl_multimin_fminimizer* s =
#        gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2rand,
#                                      model.get_num_nonlinear_params());
#
#    gsl_multimin_fminimizer_set(s, &chisq_func, x, ss);
#
#    int iter(0), status(0);
#    //
#    // Main optimization loop
#    //
#
#    do {
#        ++iter;
#        status = gsl_multimin_fminimizer_iterate(s);
#
#        if (status)
#             break;
#
#        const double size = gsl_multimin_fminimizer_size(s);
#        status = gsl_multimin_test_size(size, 1e-4);
#        //status = gsl_multimin_test_size(size, 1e-3); //changed it to test convergence 
#
#        if (status == GSL_SUCCESS)
#            std::cout << "!!! converged !!!" << std::endl;
#
#        std::cout << iter << ' ' << s->fval << ' ' << size << std::endl;
#        if (iter > 0 && iter%10 == 0) {
#            std::cout << "\\n<> solution:\\n";
#            for (size_t n = 0; n <  model.get_num_nonlinear_params(); ++n)
#                std::cout << n << " : " << s->x->data[n] << "\\n";
#            std::cout << "<>" << std::endl;
#        }
#    } while (status == GSL_CONTINUE && iter < 5000);
#
#    model.set_nonlinear_parameters(s->x->data);
#    linear::compute_chisq(s->x, 0);
#    model.set(linear::params);
#
#    gsl_vector_free(x);
#    gsl_vector_free(ss);
#    gsl_multimin_fminimizer_free(s);
#
#    //
#    // report
#    //
#
#    for (size_t i = 0; i < training_set.size(); ++i)
#      training_set[i].energy_twobody += elec_e[i] + disp_e[i];
#
#    ofstream correlation_file;
#    correlation_file.open ("correlation.dat");
#    double err_L2(0), err_wL2(0), err_Linf(0);
#    double err_L2_lo(0), err_Linf_lo(0), nlo(0);
#
#    for (size_t i = 0; i < training_set.size(); ++i) {
#
#        const double E_model = model(training_set[i].xyz) + elec_e[i] + disp_e[i] + buck_e[i];
#        const double delta = E_model - tb_ref[i];
#        if (std::abs(delta) > err_Linf)
#            err_Linf = std::abs(delta);
#
#  correlation_file <<  i+1   << "   "
#    << E_model  << "   "
#    <<  tb_ref[i]  << "    "
#          <<  delta*delta  << "    \\n" ;
#
#        err_L2 += delta*delta;
#        err_wL2 += ts_weights[i]*delta*delta;
#
#        if (training_set[i].energy_total - E_min < E_range) {
#            nlo += 1.0;
#            err_L2_lo += delta*delta;
#            if (std::abs(delta) > err_Linf_lo)
#                err_Linf_lo = std::abs(delta);
#        }
#    }
#
#    correlation_file.close();
#
#    err_L2 /= training_set.size();
#    err_wL2 /= training_set.size();
#
#    err_L2_lo /= nlo;
#
#    std::cout << "      err[L2] = " << std::sqrt(err_L2) << "    #rmsd of full ts\\n"
#              << "     err[wL2] = " << std::sqrt(err_wL2) << "   #weighted rmsd of full ts\\n"
#              << "    err[Linf] = " << err_Linf << "   #highest error in full ts\\n"
#              << "  err[L2,low] = " << std::sqrt(err_L2_lo) << "   #rmsd of low-energy ts \\n"
#              << "err[Linf,low] = " << err_Linf_lo << "   #highest error in low-energy ts "
#              << std::endl;
#
#    //
#    //  save
#    //
#
#    {
#        const char fn[] = "fit-2b.cdl";
#        std::ofstream ofs(fn);
#        std::cout << "\\n>> saving as '" << fn << "'\\n";
#        model.as_cdl(ofs);
#    }
#
#
#    return 0;
#}
#
#"""
#ff.write(a)
#ff.close()
#
#
## ## Makefile
#
## In[ ]:
#
#
#
#
## ## Evaluation code
#
## In[ ]:
#
#
#ff = open('eval-2b.cpp','w')
#a = """
##include <cmath>
##include <cassert>
##include <cstdlib>
#
##include <fstream>
##include <sstream>
##include <iomanip>
##include <iostream>
##include <stdexcept>
#
##include "mon1.h"
##include "mon2.h"
##include "training_set.h"
##include "x2b_""" + mon1 + '_' + mon2 + """_v1x.h"
##include "electrostatics.h"
##include "coulomb.h"
##include "dispersion.h"
##include "io-xyz.h"
#
#//#define GRADIENTS
#
#static std::vector<double> elec_e;
#static std::vector<double> disp_e;
#
#int main(int argc, char** argv) {
#    if (argc < 2) {
#        std::cerr << "usage: eval fit-2b.nc dimer.xyz"
#                  << std::endl;
#        return 0;
#    }
#    std::cout << std::scientific << std::setprecision(9);
#    x2b_""" + mon1 + '_' + mon2 + """::x2b_""" + mon1 + '_' + mon2 + """_v1x pot;
#    std::vector<std::string> elements;
#    std::vector<double> crd;
#
#    try {
#        ++argv;
#        --argc;
#        pot.load_netcdf(*argv);
#
#        ++argv;
#        --argc;
#        std::ifstream ifs(*argv);
#
#        if(!ifs)
#            throw std::runtime_error("could not open the XYZ file");
#
#        std::string comment;
#        kit::io::load_xyz(ifs, comment, elements, crd);
#    } catch (const std::exception& e) {
#        std::cerr << " ** Error ** : " << e.what() << std::endl;
#        return 1;
#    }
#    
#    double xyz[""" + str(3*(nat1 + nat2)) + """];
#    std::copy(crd.begin(), crd.end(), xyz);
#    
#    x::mon1 m1(xyz);
#    x::mon2 m2(xyz + 3*m1.get_realsites());
#    
#    int system_nsites = m1.get_nsites() + m2.get_nsites();
#    int * system_is_w;
#    double* system_sitecrds;
#    double* system_charge;
#    double* system_polfac;
#    double* system_pol;
#
#    system_sitecrds = new double [system_nsites*3];
#    system_charge = new double [system_nsites];
#    system_polfac = new double [system_nsites];
#    system_pol = new double [system_nsites]; //allocates memory to the pointers
#    system_is_w = new int[system_nsites];
#      
#    std::fill(system_is_w, system_is_w + m1.get_nsites(), m1.is_w);
#    std::fill(system_is_w + m1.get_nsites(), system_is_w + system_nsites, m2.is_w);
#    
#    int * is_w_a = system_is_w;
#    int * is_w_b = system_is_w + m1.get_nsites();
#    
#    std::copy(m1.get_sitecrds(), m1.get_sitecrds() + 3 * m1.get_nsites(),
#              system_sitecrds);
#    std::copy(m2.get_sitecrds(), m2.get_sitecrds() + 3 * m2.get_nsites(),
#              system_sitecrds + 3 * m1.get_nsites());
#
#    std::copy(m1.get_charges(), m1.get_charges() + m1.get_nsites(),
#              system_charge);
#    std::copy(m2.get_charges(), m2.get_charges() + m2.get_nsites(),
#              system_charge + m1.get_nsites());
#
#    std::copy(m1.get_pol(), m1.get_pol() + m1.get_nsites(),
#              system_pol);
#    std::copy(m2.get_pol(), m2.get_pol() + m2.get_nsites(),
#              system_pol + m1.get_nsites());
#
#    std::copy(m1.get_polfacs(), m1.get_polfacs() + m1.get_nsites(),
#              system_polfac);
#    std::copy(m2.get_polfacs(), m2.get_polfacs() + m2.get_nsites(),
#              system_polfac + m1.get_nsites());
#
#    excluded_set_type exclude12;
#    excluded_set_type exclude13;
#    excluded_set_type exclude14;
#
#    excluded_set_type exclude12_a;
#    excluded_set_type exclude13_a;
#    excluded_set_type exclude14_a;
#      
#    excluded_set_type exclude12_b;
#    excluded_set_type exclude13_b;
#    excluded_set_type exclude14_b;
#
#    for (auto i = m1.get_begin_12(); i != m1.get_end_12(); i++) {
#      exclude12.insert(*i);
#      exclude12_a.insert(*i);
#    }
#    for (auto i = m2.get_begin_12(); i != m2.get_end_12(); i++) {
#      std::pair<size_t,size_t> p =
#                    std::make_pair(i->first + m1.get_nsites() ,
#                                   i->second + m1.get_nsites());
#      exclude12.insert(p);
#      exclude12_b.insert(*i);
#    }
#
#    for (auto i = m1.get_begin_13(); i != m1.get_end_13(); i++) {
#      exclude13.insert(*i);
#      exclude13_a.insert(*i);
#    }
#    for (auto i = m2.get_begin_13(); i != m2.get_end_13(); i++) {
#      std::pair<size_t,size_t> p =
#                    std::make_pair(i->first + m1.get_nsites() ,
#                                   i->second + m1.get_nsites());
#      exclude13.insert(p);
#      exclude13_b.insert(*i);
#    }
#
#    for (auto i = m1.get_begin_14(); i != m1.get_end_14(); i++) {
#      exclude14.insert(*i);
#      exclude14_a.insert(*i);
#    }
#    for (auto i = m2.get_begin_14(); i != m2.get_end_14(); i++) {
#      std::pair<size_t,size_t> p =
#                    std::make_pair(i->first + m1.get_nsites() ,
#                                   i->second + m1.get_nsites());
#      exclude14.insert(p);
#      exclude14_b.insert(*i);
#    }
#
#    ttm::electrostatics m_electrostatics;
#
#    ttm::smear_ttm4x smr; 
#    smr.m_aDD_intra_12 = 0.3;
#    smr.m_aDD_intra_13 = 0.3;
#    smr.m_aDD_intra_14 = 0.055;
#
#    double ener = m_electrostatics(system_nsites, system_charge, system_polfac, system_pol,
#                                  system_sitecrds, exclude12, exclude13, exclude14, 
#                                  system_is_w, smr, 0);
#    double ener_a = m_electrostatics(m1.get_nsites(), m1.get_charges(), m1.get_polfacs(), m1.get_pol(),
#                                  m1.get_sitecrds(), exclude12_a, exclude13_a, exclude14_a, 
#                                  is_w_a, smr, 0);
#    double ener_b = m_electrostatics(m2.get_nsites(), m2.get_charges(), m2.get_polfacs(), m2.get_pol(),
#                                  m2.get_sitecrds(), exclude12_b, exclude13_b, exclude14_b, 
#                                  is_w_b, smr, 0);
#
#    // Take out electrostatic energy:
#    elec_e.push_back(ener - ener_a - ener_b);
#      
#    // Now need to take out dispersion
#    x2b_disp disp(m1.get_sitecrds(), m2.get_sitecrds(), m1.get_realsites(), m2.get_realsites());
#    ener = disp.get_dispersion();
#    disp_e.push_back(ener);
#    
#    double Epoly = pot(xyz);
#    std::cout << "IE_nograd = " << Epoly + elec_e[0] + disp_e[0] << std::endl;
#    std::cout << "E_poly2b = " << Epoly << std::endl;
#    std::cout << "E_elec2b = " << elec_e[0] << std::endl;
#    std::cout << "E_disp2b = " << disp_e[0] << std::endl;
#    
##ifdef GRADIENTS
#    const double eps = 1.0e-5;
#    double grd[""" + str(3*(nat1 + nat2)) + """];
#    std::fill(grd, grd + """ + str(3*(nat1 + nat2)) + """, 0.0);
#    Epoly = pot(xyz, grd);
#    std::cout << "E_grad = " << Epoly << std::endl;
#    for(size_t n = 0; n < """ + str(3*(nat1 + nat2)) + """; ++n){
#      const double x_orig = xyz[n];
#
#      xyz[n] = x_orig + eps;
#      const double Ep = pot(xyz) ;
#
#      xyz[n] = x_orig + 2*eps;
#      const double E2p = pot(xyz) ;
#
#      xyz[n] = x_orig - 2*eps;
#      const double E2m = pot(xyz) ;
#
#      xyz[n] = x_orig - eps;
#      const double Em = pot(xyz) ;
#
#      const double gfd = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
#      xyz[n] = x_orig;
#
#      std::cout << elements[n/3] << "   "  << "Analit: " << grd[n] << " Numerical: " << gfd
#                     << " Diff: " << std::fabs(grd[n] - gfd) << '\\n';
#    }
##endif
#
#    // Free memory
#    delete[] system_sitecrds ;
#    delete[] system_charge  ;
#    delete[] system_polfac  ;
#    delete[] system_pol ;
#    delete[] system_is_w;
#    
#    return 0;
#}
#"""
#ff.write(a)
#ff.close()
#
#
## ## Evaluation code with TTM
#
## In[ ]:
#
#
#ff = open('eval-2b-wbuck.cpp','w')
#a = """
##include <cmath>
##include <cassert>
##include <cstdlib>
#
##include <fstream>
##include <sstream>
##include <iomanip>
##include <iostream>
##include <stdexcept>
#
##include "mon1.h"
##include "mon2.h"
##include "training_set.h"
##include "x2b_""" + mon1 + '_' + mon2 + """_v1x.h"
##include "electrostatics.h"
##include "coulomb.h"
##include "dispersion.h"
##include "buckingham.h"
##include "io-xyz.h"
#
##define GRADIENTS
#
#static std::vector<double> elec_e;
#static std::vector<double> disp_e;
#static std::vector<double> buck_e;
#
#int main(int argc, char** argv) {
#    if (argc < 2) {
#        std::cerr << "usage: eval fit-2b.nc dimer.xyz"
#                  << std::endl;
#        return 0;
#    }
#    std::cout << std::scientific << std::setprecision(9);
#    x2b_""" + mon1 + '_' + mon2 + """::x2b_""" + mon1 + '_' + mon2 + """_v1x pot;
#    std::vector<std::string> elements;
#    std::vector<double> crd;
#
#    try {
#        ++argv;
#        --argc;
#        pot.load_netcdf(*argv);
#
#        ++argv;
#        --argc;
#        std::ifstream ifs(*argv);
#
#        if(!ifs)
#            throw std::runtime_error("could not open the XYZ file");
#
#        std::string comment;
#        kit::io::load_xyz(ifs, comment, elements, crd);
#    } catch (const std::exception& e) {
#        std::cerr << " ** Error ** : " << e.what() << std::endl;
#        return 1;
#    }
#    
#    double xyz[""" + str(3*(nat1 + nat2)) + """];
#    std::copy(crd.begin(), crd.end(), xyz);
#    
#    x::mon1 m1(xyz);
#    x::mon2 m2(xyz + 3*m1.get_realsites());
#    
#    int system_nsites = m1.get_nsites() + m2.get_nsites();
#    int * system_is_w;
#    double* system_sitecrds;
#    double* system_charge;
#    double* system_polfac;
#    double* system_pol;
#
#    system_sitecrds = new double [system_nsites*3];
#    system_charge = new double [system_nsites];
#    system_polfac = new double [system_nsites];
#    system_pol = new double [system_nsites]; //allocates memory to the pointers
#    system_is_w = new int[system_nsites];
#      
#    std::fill(system_is_w, system_is_w + m1.get_nsites(), m1.is_w);
#    std::fill(system_is_w + m1.get_nsites(), system_is_w + system_nsites, m2.is_w);
#    
#    int * is_w_a = system_is_w;
#    int * is_w_b = system_is_w + m1.get_nsites();
#    
#    std::copy(m1.get_sitecrds(), m1.get_sitecrds() + 3 * m1.get_nsites(),
#              system_sitecrds);
#    std::copy(m2.get_sitecrds(), m2.get_sitecrds() + 3 * m2.get_nsites(),
#              system_sitecrds + 3 * m1.get_nsites());
#
#    std::copy(m1.get_charges(), m1.get_charges() + m1.get_nsites(),
#              system_charge);
#    std::copy(m2.get_charges(), m2.get_charges() + m2.get_nsites(),
#              system_charge + m1.get_nsites());
#
#    std::copy(m1.get_pol(), m1.get_pol() + m1.get_nsites(),
#              system_pol);
#    std::copy(m2.get_pol(), m2.get_pol() + m2.get_nsites(),
#              system_pol + m1.get_nsites());
#
#    std::copy(m1.get_polfacs(), m1.get_polfacs() + m1.get_nsites(),
#              system_polfac);
#    std::copy(m2.get_polfacs(), m2.get_polfacs() + m2.get_nsites(),
#              system_polfac + m1.get_nsites());
#
#    excluded_set_type exclude12;
#    excluded_set_type exclude13;
#    excluded_set_type exclude14;
#
#    excluded_set_type exclude12_a;
#    excluded_set_type exclude13_a;
#    excluded_set_type exclude14_a;
#      
#    excluded_set_type exclude12_b;
#    excluded_set_type exclude13_b;
#    excluded_set_type exclude14_b;
#
#    for (auto i = m1.get_begin_12(); i != m1.get_end_12(); i++) {
#      exclude12.insert(*i);
#      exclude12_a.insert(*i);
#    }
#    for (auto i = m2.get_begin_12(); i != m2.get_end_12(); i++) {
#      std::pair<size_t,size_t> p =
#                    std::make_pair(i->first + m1.get_nsites() ,
#                                   i->second + m1.get_nsites());
#      exclude12.insert(p);
#      exclude12_b.insert(*i);
#    }
#
#    for (auto i = m1.get_begin_13(); i != m1.get_end_13(); i++) {
#      exclude13.insert(*i);
#      exclude13_a.insert(*i);
#    }
#    for (auto i = m2.get_begin_13(); i != m2.get_end_13(); i++) {
#      std::pair<size_t,size_t> p =
#                    std::make_pair(i->first + m1.get_nsites() ,
#                                   i->second + m1.get_nsites());
#      exclude13.insert(p);
#      exclude13_b.insert(*i);
#    }
#
#    for (auto i = m1.get_begin_14(); i != m1.get_end_14(); i++) {
#      exclude14.insert(*i);
#      exclude14_a.insert(*i);
#    }
#    for (auto i = m2.get_begin_14(); i != m2.get_end_14(); i++) {
#      std::pair<size_t,size_t> p =
#                    std::make_pair(i->first + m1.get_nsites() ,
#                                   i->second + m1.get_nsites());
#      exclude14.insert(p);
#      exclude14_b.insert(*i);
#    }
#
#    ttm::electrostatics m_electrostatics;
#
#    ttm::smear_ttm4x smr; 
#    smr.m_aDD_intra_12 = 0.3;
#    smr.m_aDD_intra_13 = 0.3;
#    smr.m_aDD_intra_14 = 0.055;
#
#    double ener = m_electrostatics(system_nsites, system_charge, system_polfac, system_pol,
#                                  system_sitecrds, exclude12, exclude13, exclude14, 
#                                  system_is_w, smr, 0);
#    double ener_a = m_electrostatics(m1.get_nsites(), m1.get_charges(), m1.get_polfacs(), m1.get_pol(),
#                                  m1.get_sitecrds(), exclude12_a, exclude13_a, exclude14_a, 
#                                  is_w_a, smr, 0);
#    double ener_b = m_electrostatics(m2.get_nsites(), m2.get_charges(), m2.get_polfacs(), m2.get_pol(),
#                                  m2.get_sitecrds(), exclude12_b, exclude13_b, exclude14_b, 
#                                  is_w_b, smr, 0);
#
#    // Take out electrostatic energy:
#    elec_e.push_back(ener - ener_a - ener_b);
#      
#    // Now need to take out dispersion
#    x2b_disp disp(m1.get_sitecrds(), m2.get_sitecrds(), m1.get_realsites(), m2.get_realsites());
#    ener = disp.get_dispersion();
#    disp_e.push_back(ener);
#    
#    // Now need to take out buckingham
#    x2b_buck buck(m1.get_sitecrds(), m2.get_sitecrds(), m1.get_realsites(), m2.get_realsites());
#    ener = buck.get_buckingham();
#    buck_e.push_back(ener);
#    
#    double Epoly = pot(xyz);
#    std::cout << "IE_nograd = " << Epoly + elec_e[0] + disp_e[0] << std::endl;
#    std::cout << "E_poly2b = " << Epoly << std::endl;
#    std::cout << "E_elec2b = " << elec_e[0] << std::endl;
#    std::cout << "E_disp2b = " << disp_e[0] << std::endl;
#    std::cout << "E_buck2b = " << buck_e[0] << std::endl;
#    
##ifdef GRADIENTS
#    const double eps = 1.0e-5;
#    double grd[""" + str(3*(nat1 + nat2)) + """];
#    std::fill(grd, grd + """ + str(3*(nat1 + nat2)) + """, 0.0);
#    Epoly = pot(xyz, grd);
#    std::cout << "E_grad = " << Epoly << std::endl;
#    for(size_t n = 0; n < """ + str(3*(nat1 + nat2)) + """; ++n){
#      const double x_orig = xyz[n];
#
#      xyz[n] = x_orig + eps;
#      const double Ep = pot(xyz) ;
#
#      xyz[n] = x_orig + 2*eps;
#      const double E2p = pot(xyz) ;
#
#      xyz[n] = x_orig - 2*eps;
#      const double E2m = pot(xyz) ;
#
#      xyz[n] = x_orig - eps;
#      const double Em = pot(xyz) ;
#
#      const double gfd = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
#      xyz[n] = x_orig;
#
#      std::cout << elements[n/3] << "   "  << "Analit: " << grd[n] << " Numerical: " << gfd
#                     << " Diff: " << std::fabs(grd[n] - gfd) << '\\n';
#    }
##endif
#
#    // Free memory
#    delete[] system_sitecrds ;
#    delete[] system_charge  ;
#    delete[] system_polfac  ;
#    delete[] system_pol ;
#    delete[] system_is_w;
#    
#    return 0;
#}
#"""
#ff.write(a)
#ff.close()
#
#
## ## X2B.h for eval    
#
## In[ ]:
#
#
#hname = "x2b_" + mon1 + "_" + mon2 + "_v1x.h"
#polyhname = "poly_2b_" + mon1 + "_" + mon2 + "_v1x.h"
#defname = "X2B_" + mon1 + "_" + mon2 + "_V1X_H"
#ff = open(hname,'w')
#ff.write('#ifndef ' + defname + '\n')
#ff.write('#define ' + defname + '\n \n')
#ff.write('#include "' + polyhname + '" \n')
#
#a = """
##include <iostream>
#////////////////////////////////////////////////////////////////////////////////
#
#namespace x2b_""" + mon1 + "_" + mon2 + """ {
#
#//----------------------------------------------------------------------------//
#
#"""
#ff.write(a)
#ff.write('struct x2b_' + mon1 + "_" + mon2 + "_v1x { \n")
#a = """
#    typedef mb_system::poly_model poly_type;
#
#
#    static std::string name();
#    void load_netcdf(const char*);
#
#    // returns 2B contribution only
#    // XYZ is for the real sites
#"""
#ff.write(a)
#ff.write('    double operator()(const double xyz[' + str(3*(nat1 + nat2)) + '], double grd[' + str(3*(nat1 + nat2)) + ']) const; \n')
#ff.write('    double operator()(const double xyz[' + str(3*(nat1 + nat2)) + ']) const; \n')
#a = """
#    double eval(const double* mon1, const double* mon2) const;
#    double eval(const double* mon1, const double* mon2, double* g1, double* g2) const;
#    
#private:
#
#"""
#ff.write(a)
#for nl in nlparam_all:
#    ff.write("    double m_" + nl + ";\n")
#a = """
#protected:
#    double m_r2i = """ + str(r2i) + """;
#    double m_r2f = """ + str(r2o) + """;
#
#    double f_switch(const double&, double&) const; // At1_a -- At1_b separation
#
#private:
#    double m_poly[poly_type::size];
#};
#
#//----------------------------------------------------------------------------//
#
#} // namespace x2b_""" + mon1 + "_" + mon2 + """
#
#////////////////////////////////////////////////////////////////////////////////
#
##endif 
#"""
#ff.write(a)
#ff.close()
#
#
## ## X2B.cpp for eval    
#
## In[ ]:
#
#
#cppname = "x2b_" + mon1 + "_" + mon2 + "_v1x.cpp"
#ff = open(cppname,'w')
#a = """
##include <cmath>
##include <cassert>
##include <cstdlib>
##include <iomanip>
#
##include <netcdf.h>
#
#"""
#hname = "x2b_" + mon1 + "_" + mon2 + "_v1x.h"
#ff.write(a)
#ff.write('#include "' + hname + '" \n \n')
#a = """
#////////////////////////////////////////////////////////////////////////////////
#
#namespace {
#
#//----------------------------------------------------------------------------//
#
#void error(int kode) {
#
#    std::cerr << " ** Fatal Error in x2b_""" + mon1 + "_" + mon2 + """_v1x::load_netcdf() ** :" 
#              << nc_strerror(kode) << std::endl;
#    std::exit(EXIT_FAILURE);
#}
#
#//----------------------------------------------------------------------------//
#
#struct variable {
#    double v_exp0(const double& r0, const double& k,
#                 const double * p1, const double * p2 );
#                 
#    double v_exp(const double& k,
#                 const double * p1, const double * p2 );
#
#    double v_coul0(const double& r0, const double& k,
#                  const double * p1, const double * p2 );
#                  
#    double v_coul(const double& k,
#                  const double * p1, const double * p2 );
#
#    double v_gau0(const double& r0, const double& k,
#                 const double * p1, const double * p2 );
#                  
#    void grads(const double& gg, double * grd1, double * grd2,
#               const double * p1, const double * p2);
#
#    double g[3]; // diff(value, p1 - p2)
#};
#
#//----------------------------------------------------------------------------//
#
#double variable::v_gau0(const double& r0, const double& k,
#                       const double * p1, const double * p2)
#{
#    g[0] = p1[0] - p2[0];
#    g[1] = p1[1] - p2[1];
#    g[2] = p1[2] - p2[2];
#
#    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
#
#    const double exp1 = std::exp(-k*(r0 - r)*(r0 - r));
#    const double gg = 2*k*(r0 - r)*exp1/r;
#
#    g[0] *= gg;
#    g[1] *= gg;
#    g[2] *= gg;
#
#    return exp1;
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#double variable::v_exp0(const double& r0, const double& k,
#                       const double * p1, const double * p2)
#{
#    g[0] = p1[0] - p2[0];
#    g[1] = p1[1] - p2[1];
#    g[2] = p1[2] - p2[2];
#
#    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
#
#    const double exp1 = std::exp(k*(r0 - r));
#    const double gg = - k*exp1/r;
#
#    g[0] *= gg;
#    g[1] *= gg;
#    g[2] *= gg;
#
#    return exp1;
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#double variable::v_exp(const double& k,
#                       const double * p1, const double * p2)
#{
#    g[0] = p1[0] - p2[0];
#    g[1] = p1[1] - p2[1];
#    g[2] = p1[2] - p2[2];
#
#    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
#
#    const double exp1 = std::exp(k*(- r));
#    const double gg = - k*exp1/r;
#
#    g[0] *= gg;
#    g[1] *= gg;
#    g[2] *= gg;
#
#    return exp1;
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#
#double variable::v_coul(const double& k,
#                        const double * p1, const double * p2)
#{
#    g[0] = p1[0] - p2[0];
#    g[1] = p1[1] - p2[1];
#    g[2] = p1[2] - p2[2];
#
#    const double rsq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
#    const double r = std::sqrt(rsq);
#
#    const double exp1 = std::exp(k*(-r));
#    const double rinv = 1.0/r;
#    const double val = exp1*rinv;
#
#    const double gg = - (k + rinv)*val*rinv;
#
#    g[0] *= gg;
#    g[1] *= gg;
#    g[2] *= gg;
#
#    return val;
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#double variable::v_coul0(const double& r0, const double& k,
#                        const double * p1, const double * p2)
#{
#    g[0] = p1[0] - p2[0];
#    g[1] = p1[1] - p2[1];
#    g[2] = p1[2] - p2[2];
#
#    const double rsq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
#    const double r = std::sqrt(rsq);
#
#    const double exp1 = std::exp(k*(r0 - r));
#    const double rinv = 1.0/r;
#    const double val = exp1*rinv;
#
#    const double gg = - (k + rinv)*val*rinv;
#
#    g[0] *= gg;
#    g[1] *= gg;
#    g[2] *= gg;
#
#    return val;
#}
#
#//----------------------------------------------------------------------------//
#
#void variable::grads(const double& gg, double * grd1, double * grd2, 
#                     const double * p1, const double * p2) {
#    for (size_t i = 0; i < 3 ; i++) {
#        double d = gg*g[i];
#        grd1[i] += d;
#        grd2[i] -= d;
#    }
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#struct monomer {
#    double oh1[3];
#    double oh2[3];
#
#    void setup(const double* ohh,
#               const double& in_plane_g, const double& out_of_plane_g,
#               double x1[3], double x2[3]);
#
#    void grads(const double* g1, const double* g2,
#               const double& in_plane_g, const double& out_of_plane_g,
#               double* grd) const;
#};
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#void monomer::setup(const double* ohh,
#                    const double& in_plane_g, const double& out_of_plane_g,
#                    double* x1, double* x2)
#{
#    for (int i = 0; i < 3; ++i) {
#        oh1[i] = ohh[i + 3] - ohh[i];
#        oh2[i] = ohh[i + 6] - ohh[i];
#    }
#
#    const double v[3] = {
#        oh1[1]*oh2[2] - oh1[2]*oh2[1],
#        oh1[2]*oh2[0] - oh1[0]*oh2[2],
#        oh1[0]*oh2[1] - oh1[1]*oh2[0]
#    };
#
#    for (int i = 0; i < 3; ++i) {
#        const double in_plane = ohh[i] + 0.5*in_plane_g*(oh1[i] + oh2[i]);
#        const double out_of_plane = out_of_plane_g*v[i];
#
#        x1[i] = in_plane + out_of_plane;
#        x2[i] = in_plane - out_of_plane;
#    }
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#void monomer::grads(const double* g1, const double* g2,
#                    const double& in_plane_g, const double& out_of_plane_g,
#                    double* grd) const
#{
#    const double gm[3] = {
#        g1[0] - g2[0],
#        g1[1] - g2[1],
#        g1[2] - g2[2]
#    };
#
#    const double t1[3] = {
#        oh2[1]*gm[2] - oh2[2]*gm[1],
#        oh2[2]*gm[0] - oh2[0]*gm[2],
#        oh2[0]*gm[1] - oh2[1]*gm[0]
#    };
#
#    const double t2[3] = {
#        oh1[1]*gm[2] - oh1[2]*gm[1],
#        oh1[2]*gm[0] - oh1[0]*gm[2],
#        oh1[0]*gm[1] - oh1[1]*gm[0]
#    };
#
#    for (int i = 0; i < 3; ++i) {
#        const double gsum = g1[i] + g2[i];
#        const double in_plane = 0.5*in_plane_g*gsum;
#
#        const double gh1 = in_plane + out_of_plane_g*t1[i];
#        const double gh2 = in_plane - out_of_plane_g*t2[i];
#
#        grd[i + 0] += gsum - (gh1 + gh2); // O
#        grd[i + 3] += gh1; // H1
#        grd[i + 6] += gh2; // H2
#    }
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#//struct vsites {
#//    //void TwoParticleAverageSite() {}
#//    //void ThreeParticleAverageSite() {}
#//    void OutOfPlaneSite(const double& w12, const double& w13,
#//                        const double& wcross, const double x1[3],
#//                        const double y1[3], const double y2[3],
#//                        double vs[3]);
#//    //void LocalCoordinatesSite{}
#//};
#//
#//void vsites::OutOfPlaneSite(const double& w12,
#//                            const double& w13,
#//                            const double& wcross,
#//                            const double x1[3],
#//                            const double y1[3],
#//                            const double y2[3],
#//                            double vs[3]) {
#//    double r12[3], r13[3];
#//
#//    for (int i = 0; i < 3; ++i) {
#//        r12[i] = y1[i] - x1[i];
#//        r13[i] = y2[i] - x1[i];
#//    }
#//                            
#//    double rc[3];
#//    rc[0] = r12[1]*r13[2] - r12[2]*r13[1];
#//    rc[1] = r12[2]*r13[0] - r12[0]*r13[2];
#//    rc[2] = r12[0]*r13[1] - r12[1]*r13[0];
#//    
#//    vs[0] = x1[0] + w12 * r12[0] + w13 * r13[0] + wcross * rc[0];
#//    vs[1] = x1[1] + w12 * r12[1] + w13 * r13[1] + wcross * rc[1];
#//    vs[2] = x1[2] + w12 * r12[2] + w13 * r13[2] + wcross * rc[2];
#//}
#
#} // namespace
#
#////////////////////////////////////////////////////////////////////////////////
#
#namespace x2b_""" + mon1 + "_" + mon2 + """ {
#
#//----------------------------------------------------------------------------//
#
#
#std::string x2b_""" + mon1 + "_" + mon2 + """_v1x::name() {
#    return "x2b_""" + mon1 + "_" + mon2 + """_v1x";
#}
#
#"""
#ff.write(a)
#a = """
#//----------------------------------------------------------------------------//
#
#void x2b_""" + mon1 + "_" + mon2 + """_v1x::load_netcdf(const char* fn)
#{
#    assert(fn);
#
#    int rc, ncid;
#    if ((rc = nc_open(fn, NC_NOWRITE, &ncid)))
#        error(rc);
#
##   define RETRIEVE(name) \\
#    if ((rc = nc_get_att_double(ncid, NC_GLOBAL, #name , &m_##name))) \\
#        error(rc);
#"""
#ff.write(a)
#for nl in nlparam_all:
#    ff.write("    RETRIEVE(" + nl + ")\n")
#a = """
#
#    RETRIEVE(r2i)
#    RETRIEVE(r2f)
#
##   undef RETRIEVE
#
#    int varid;
#
#    if ((rc = nc_inq_varid(ncid, "poly", &varid)))
#        error(rc);
#
#    for (size_t n = 0; n < poly_type::size; ++n) {
#        if ((rc = nc_get_var1_double(ncid, varid, &n, m_poly + n)))
#            error(rc);
#    }
#
#    if ((rc = nc_close(ncid)))
#        error(rc);
#}
#
#//----------------------------------------------------------------------------//
#
#double x2b_""" + mon1 + "_" + mon2 + """_v1x::f_switch(const double& r, double& g) const
#{
#    if (r > m_r2f) {
#        g = 0.0;
#        return 0.0;
#    } else if (r > m_r2i) {
#        const double t1 = M_PI/(m_r2f - m_r2i);
#        const double x = (r - m_r2i)*t1;
#        g = - std::sin(x)*t1/2.0;
#        return (1.0 + std::cos(x))/2.0;
#    } else {
#        g = 0.0;
#        return 1.0;
#    }
#}
#
#//----------------------------------------------------------------------------//
#
#double x2b_""" + mon1 + "_" + mon2 + """_v1x::eval(const double* mon1, const double* mon2 ) const
#{
#    // the switch
#
#    // ##DEFINE HERE## right now it assumes 1st atom of each monomer
#    const double d12[3] = {mon1[0] -  mon2[0],
#                           mon1[1] -  mon2[1],
#                           mon1[2] -  mon2[2]};
#
#    const double r12sq = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
#    const double r12 = std::sqrt(r12sq);
#
#    if (r12 > m_r2f)
#        return 0.0;
#
#    double xcrd[""" + str(3*(nat1 + nat2)) + """]; // coordinates of real sites ONLY
#
#    std::copy(mon1, mon1 + """ + str(3*(nat1)) + """, xcrd);
#    std::copy(mon2, mon2 + """ + str(3*(nat2)) + """, xcrd + """ + str(3*(nat1)) + """);
#    
#    double v[""" + str(nvars) + """];
#    
#    double sw = 0.0;
#    double gsw = 0.0;
#    
#"""
#ff.write(a)
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_a),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_a[i+1])):
#        if not types_a[i] in vsites:
#            ff.write('    const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xcrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_a),2):
#    n = 1
#    for j in range(int(types_a[i+1])):
#        if types_a[i] in vsites:
#            ff.write('    double ' + types_a[i] + '_' + str(n) + '_a[3]' + ';\n')
#            n = n + 1
#ff.write('\n')
#
## loops over each type of atom in the input
#for i in range(0,len(types_b),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_b[i+1])):
#        if not types_b[i] in vsites:
#            ff.write('    const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xcrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_b),2):
#    n = 1
#    for j in range(int(types_b[i+1])):
#        if types_b[i] in vsites:
#            ff.write('    double ' + types_b[i] + '_' + str(n) + '_b[3]' + ';\n')
#            n = n + 1
#ff.write('\n')
#
#if use_lonepairs[0] == 0 and use_lonepairs[1] == 0:
#    a = """    
#    variable vr[""" + str(nvars) + """];
#    
#"""
#elif use_lonepairs[0] != 0:
#    a = """
#//    vsites virt;
#    double w12 =     -9.721486914088159e-02;  //from MBpol
#    double w13 =     -9.721486914088159e-02;
#    double wcross =   9.859272078406150e-02;
#
#    monomer m;
#    
#    m.setup(""" + types_a[0] + '_1_a' + """, w12, wcross,
#             """ + types_a[4] + '_1_a' + """, """ + types_a[4] + '_2_a' + """);
#                    
#    variable vr[""" + str(nvars) + """];
#    
#"""
#elif use_lonepairs[1] != 0:
#    a = """
#//    vsites virt;
#    double w12 =     -9.721486914088159e-02;  //from MBpol
#    double w13 =     -9.721486914088159e-02;
#    double wcross =   9.859272078406150e-02;
#
#    monomer m;
#    
#    m.setup(""" + types_b[0] + '_1_b' + """, w12, wcross,
#             """ + types_b[4] + '_1_b' + """, """ + types_b[4] + '_2_b' + """);
#                    
#    variable vr[""" + str(nvars) + """];
#    
#"""
#ff.write(a)
#
#nv = 0
## Intramolecular distances:
#for index, variable in enumerate(variables):
#    if variable[4].startswith("x-intra-"):
#        atom1 = variable[0][0]
#        atom2 = variable[2][0]
#
#        try:
#            atom1_index = int(variable[0][1:])
#        except ValueError:
#            atom1_index = 1
#
#        try:
#            atom2_index = int(variable[2][1:])
#        except ValueError:
#            atom2_index = 1
#
#        atom1_fragment = variable[1]
#        atom2_fragment = variable[3]
#
#        atom1_name = "{}_{}_{}".format(atom1, atom1_index, atom1_fragment)
#        atom2_name = "{}_{}_{}".format(atom2, atom2_index, atom2_fragment)
#
#        sorted_atoms = "".join(sorted([atom1, atom2]))
#
#        if var_intra == 'exp0' or var_intra == 'coul0' or var_intra == 'gau0':
#            arguments = '(m_d_intra_' + sorted_atoms + ', m_k_intra_' + sorted_atoms
#        else:
#            arguments = '(m_k_intra_' + sorted_atoms
#
#        ff.write('    v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var_intra + arguments + ', ' + atom1_name + ', ' + atom2_name + ');\n')
#
#        nv += 1
#    else:
#        # inter-molecular variables
#        atom1 = variable[0][0]
#        atom2 = variable[2][0]
#
#        try:
#            atom1_index = int(variable[0][1:])
#        except ValueError:
#            atom1_index = 1
#
#        try:
#            atom2_index = int(variable[2][1:])
#        except ValueError:
#            atom2_index = 1
#
#        atom1_fragment = variable[1]
#        atom2_fragment = variable[3]
#
#        atom1_name = "{}_{}_{}".format(atom1, atom1_index, atom1_fragment)
#        atom2_name = "{}_{}_{}".format(atom2, atom2_index, atom2_fragment)
#
#        if atom1 not in vsites and atom2 not in vsites:
#            var = var_inter
#        else:
#            var = var_lp
#
#        sorted_atoms = "".join(sorted([atom1, atom2]))
#
#        if var == 'exp0' or var == 'coul0' or var == 'gau0':
#            arguments = '(m_d_' + sorted_atoms + ', m_k_' + sorted_atoms
#        else:
#            arguments = '(m_k_' + sorted_atoms
#
#        ff.write('    v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var + arguments + ', ' + atom1_name + ', ' + atom2_name + ');\n')
#
#        nv += 1
# 
#a = """     
#    
#    sw = f_switch(r12, gsw);
#    
#    const double E_poly = mb_system::poly_model::eval(m_poly, v);
#    
#    return sw*E_poly;
#}
#
#double x2b_""" + mon1 + "_" + mon2 + """_v1x::eval(const double* mon1, const double* mon2, 
#                double * grd1, double * grd2) const
#{
#
#    // ##DEFINE HERE## right now it assumes 1st atom of each monomer
#    const double d12[3] = {mon1[0] -  mon2[0],
#                           mon1[1] -  mon2[1],
#                           mon1[2] -  mon2[2]};
#
#    const double r12sq = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
#    const double r12 = std::sqrt(r12sq);
#
#    if (r12 > m_r2f)
#        return 0.0;
#
#    double xcrd[""" + str(3*(nat1 + nat2)) + """]; // coordinates of real sites ONLY
#
#    std::copy(mon1, mon1 + """ + str(3*(nat1)) + """, xcrd);
#    std::copy(mon2, mon2 + """ + str(3*(nat2)) + """, xcrd + """ + str(3*(nat1)) + """);
#    
#    double v[""" + str(nvars) + """];
#    
#    double sw = 0.0;
#    double gsw = 0.0;
#    
#"""
#ff.write(a)
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_a),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_a[i+1])):
#        if not types_a[i] in vsites:
#            ff.write('    const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xcrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_a),2):
#    n = 1
#    for j in range(int(types_a[i+1])):
#        if types_a[i] in vsites:
#            ff.write('    double ' + types_a[i] + '_' + str(n) + '_a[3]' + ';\n')
#            n = n + 1
#ff.write('\n')
#
## loops over each type of atom in the input
#for i in range(0,len(types_b),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_b[i+1])):
#        if not types_b[i] in vsites:
#            ff.write('    const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xcrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_b),2):
#    n = 1
#    for j in range(int(types_b[i+1])):
#        if types_b[i] in vsites:
#            ff.write('    double ' + types_b[i] + '_' + str(n) + '_b[3]' + ';\n')
#            n = n + 1
#ff.write('\n')
#
#if use_lonepairs[0] == 0 and use_lonepairs[1] == 0:
#    a = """    
#    variable vr[""" + str(nvars) + """];
#    
#"""
#elif use_lonepairs[0] != 0:
#    a = """
#    //vsites virt;
#    double w12 =     -9.721486914088159e-02;  //from MBpol
#    double w13 =     -9.721486914088159e-02;
#    double wcross =   9.859272078406150e-02;
#
#    monomer m;
#    
#    m.setup(""" + types_a[0] + '_1_a' + """, w12, wcross, 
#             """ + types_a[4] + '_1_a' + """, """ + types_a[4] + '_2_a' + """);
#                    
#    variable vr[""" + str(nvars) + """];
#    
#"""
#elif use_lonepairs[1] != 0:
#    a = """
#    //vsites virt;
#    double w12 =     -9.721486914088159e-02;  //from MBpol
#    double w13 =     -9.721486914088159e-02;
#    double wcross =   9.859272078406150e-02;
#
#    monomer m;
#    
#    m.setup(""" + types_b[0] + '_1_b' + """, w12, wcross, 
#             """ + types_b[4] + '_1_b' + """, """ + types_b[4] + '_2_b' + """);
#                    
#    variable vr[""" + str(nvars) + """];
#    
#"""
#ff.write(a)
#
#nv = 0
## Intramolecular distances:
#for index, variable in enumerate(variables):
#    if variable[4].startswith("x-intra-"):
#        atom1 = variable[0][0]
#        atom2 = variable[2][0]
#
#        try:
#            atom1_index = int(variable[0][1:])
#        except ValueError:
#            atom1_index = 1
#
#        try:
#            atom2_index = int(variable[2][1:])
#        except ValueError:
#            atom2_index = 1
#
#        atom1_fragment = variable[1]
#        atom2_fragment = variable[3]
#
#        atom1_name = "{}_{}_{}".format(atom1, atom1_index, atom1_fragment)
#        atom2_name = "{}_{}_{}".format(atom2, atom2_index, atom2_fragment)
#
#        sorted_atoms = "".join(sorted([atom1, atom2]))
#
#        if var_intra == 'exp0' or var_intra == 'coul0' or var_intra == 'gau0':
#            arguments = '(m_d_intra_' + sorted_atoms + ', m_k_intra_' + sorted_atoms
#        else:
#            arguments = '(m_k_intra_' + sorted_atoms
#
#        ff.write('    v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var_intra + arguments + ', ' + atom1_name + ', ' + atom2_name + ');\n')
#
#        nv += 1
#    else:
#        # inter-molecular variables
#        atom1 = variable[0][0]
#        atom2 = variable[2][0]
#
#        try:
#            atom1_index = int(variable[0][1:])
#        except ValueError:
#            atom1_index = 1
#
#        try:
#            atom2_index = int(variable[2][1:])
#        except ValueError:
#            atom2_index = 1
#
#        atom1_fragment = variable[1]
#        atom2_fragment = variable[3]
#
#        atom1_name = "{}_{}_{}".format(atom1, atom1_index, atom1_fragment)
#        atom2_name = "{}_{}_{}".format(atom2, atom2_index, atom2_fragment)
#
#        if atom1 not in vsites and atom2 not in vsites:
#            var = var_inter
#        else:
#            var = var_lp
#
#        sorted_atoms = "".join(sorted([atom1, atom2]))
#
#        if var == 'exp0' or var == 'coul0' or var == 'gau0':
#            arguments = '(m_d_' + sorted_atoms + ', m_k_' + sorted_atoms
#        else:
#            arguments = '(m_k_' + sorted_atoms
#
#        ff.write('    v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var + arguments + ', ' + atom1_name + ', ' + atom2_name + ');\n')
#
#        nv += 1
#
#a = """     
#    
#    double g[""" + str(nvars) + """];
#    
#    const double E_poly = mb_system::poly_model::eval(m_poly, v, g);
#    
#    double xgrd[""" + str(3*(len(atom_list_a) + len(atom_list_b))) + """];
#    std::fill(xgrd, xgrd + """ + str(3*(len(atom_list_a) + len(atom_list_b))) + """, 0.0);
#
#""" 
#ff.write(a)   
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_a),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_a[i+1])):
#        if not types_a[i] in vsites:
#            ff.write('    double* ' + types_a[i] + '_' + str(n) + '_a_g' + '= xgrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
## loops over each type of atom in the input
#for i in range(0,len(types_b),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_b[i+1])):
#        if not types_b[i] in vsites:
#            ff.write('    double* ' + types_b[i] + '_' + str(n) + '_b_g' + '= xgrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_a),2):
#    n = 1
#    for j in range(int(types_a[i+1])):
#        if types_a[i] in vsites:
#            ff.write('    double* ' + types_a[i] + '_' + str(n) + '_a_g' + '= xgrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#ff.write('\n')
#
#for i in range(0,len(types_b),2):
#    n = 1
#    for j in range(int(types_b[i+1])):
#        if types_b[i] in vsites:
#            ff.write('    double* ' + types_b[i] + '_' + str(n) + '_b_g' + '= xgrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#ff.write('\n')
#
#nv = 0
## Intramolecular distances:
#for index, variable in enumerate(variables):
#    atom1 = variable[0][0]
#    atom2 = variable[2][0]
#
#    try:
#        atom1_index = int(variable[0][1:])
#    except ValueError:
#        atom1_index = 1
#
#    try:
#        atom2_index = int(variable[2][1:])
#    except ValueError:
#        atom2_index = 1
#
#    atom1_fragment = variable[1]
#    atom2_fragment = variable[3]
#
#    atom1_name = "{}_{}_{}".format(atom1, atom1_index, atom1_fragment)
#    atom2_name = "{}_{}_{}".format(atom2, atom2_index, atom2_fragment)
#
#    ff.write('    vr[' + str(nv) + '].grads(g[' + str(nv) + '], ' + atom1_name + '_g, ' + atom2_name + '_g, ' + atom1_name + ', ' + atom2_name + ');\n')
#
#    nv += 1
#            
#a = """
#
#    // ##DEFINE HERE## the redistribution of the gradients
#    
#"""
#ff.write(a)
#a = ""
#if use_lonepairs[0] != 0:
#    a = """
#    m.grads(""" + types_a[4] + '_1_a_g, ' + types_a[4] + """_2_a_g, 
#             w12, wcross, """ + types_a[0] + """_1_a_g);
#    """
#elif use_lonepairs[1] != 0:
#    a = """
#    m.grads(""" + types_b[4] + '_1_b_g, ' + types_b[4] + """_2_b_g, 
#             w12, wcross, """ + types_b[0] + """_1_b_g);
#"""
#ff.write(a)    
#
#a = """
#    
#    // the switch
#    
#    sw = f_switch(r12, gsw);
#    
#    for (int i = 0; i < """ + str(3*nat1) + """; ++i) {
#        grd1[i] += sw*xgrd[i];
#    }
#
#    for (int i = 0; i < """ + str(3*nat2) + """; ++i) {
#        grd2[i] += sw*xgrd[i + """ + str(3*nat1) + """];
#    }
#
#    // gradient of the switch
#
#    gsw *= E_poly/r12;
#    for (int i = 0; i < 3; ++i) {
#        const double d = gsw*d12[i];
#        grd1[i] += d;
#        grd2[i] -= d;
#    }
#    return sw*E_poly;
#}
#
#double x2b_""" + mon1 + "_" + mon2 + """_v1x::operator()(const double crd[""" + str(3*(nat1 + nat2)) + """]) const
#{
#    const double E_poly = eval(crd, crd + """ + str(3*(nat1)) + """);
#
#    return E_poly;
#}
#
#double x2b_""" + mon1 + "_" + mon2 + """_v1x::operator()(const double crd[""" + str(3*(nat1 + nat2)) + """],
#                            double grd[""" + str(3*(nat1 + nat2)) + """]) const
#{
#    const double E_poly = eval(crd, crd + """ + str(3*(nat1)) + """, grd, grd + """ + str(3*(nat1)) + """);
#
#    return E_poly;
#}
#
#} // namespace x2b_""" + mon1 + "_" + mon2 + """
#
#////////////////////////////////////////////////////////////////////////////////
#"""
#ff.write(a)
#ff.close()
#
#
## ## X2B.h for software
#
## In[ ]:
#
#
#hname = "x2b_" + mon1 + "_" + mon2 + "_deg" + str(degree) + "_v1x.h"
#polyhname = "poly_2b_" + mon1 + "_" + mon2 + "_deg" + str(degree) + "_v1x.h"
#defname = "X2B_" + mon1 + "_" + mon2 + "_DEG" + str(degree) + "_V1X_H"
#ff = open(hname,'w')
#ff.write('#ifndef ' + defname + '\n')
#ff.write('#define ' + defname + '\n \n')
#ff.write('#include "' + polyhname + '" \n')
#
#a = """
##include <iostream>
##include <cmath>
##include <cassert>
##include <cstdlib>
##include <iomanip>
##include <string>
##include <vector>
#
#////////////////////////////////////////////////////////////////////////////////
#
#namespace x2b_""" + mon1 + "_" + mon2 + "_deg" + str(degree) + """ {
#
#//----------------------------------------------------------------------------//
#
#"""
#ff.write(a)
#ff.write('struct x2b_' + mon1 + "_" + mon2 + "_v1x { \n")
#a = """
#    x2b_""" + mon1 + "_" + mon2 + """_v1x() {};
#    x2b_""" + mon1 + "_" + mon2 + """_v1x(const std::string mon1, const std::string mon2);
#
#    ~x2b_""" + mon1 + "_" + mon2 + """_v1x() {};
#
#    typedef poly_2b_""" + mon1 + "_" + mon2 + "_deg" + str(degree) + """_v1x polynomial;
#
#    // returns 2B contribution only
#    // XYZ is for the real sites
#"""
#ff.write(a)
#a = """
#    double eval(const double* xyz1, const double* xyz2, const size_t ndim) const;
#    double eval(const double* xyz1, const double* xyz2, double* grad1, double* grad2, const size_t ndim) const;
#    
#private:
#
#"""
#ff.write(a)
#for nl in nlparam_all:
#    ff.write("    double m_" + nl + ";\n")
#a = """
#protected:
#    double m_r2i = """ + str(r2i) + """;
#    double m_r2f = """ + str(r2o) + """;
#
#    double f_switch(const double&, double&) const; // At1_a -- At1_b separation
#
#private:
#    std::vector<double> coefficients;
#};
#
#//----------------------------------------------------------------------------//
#
#} // namespace x2b_""" + mon1 + "_" + mon2 + "_deg" + str(degree) + """
#
#////////////////////////////////////////////////////////////////////////////////
#
##endif 
#"""
#ff.write(a)
#ff.close()
#
#
## ## X2B.cpp for software
#
## In[ ]:
#
#
#cppname = "x2b_" + mon1 + "_" + mon2 + "_deg" + str(degree) + "_v1x.cpp"
#ff = open(cppname,'w')
#a = """
#
#"""
#hname = "x2b_" + mon1 + "_" + mon2 + "_deg" + str(degree) + "_v1x.h"
#ff.write(a)
#ff.write('#include "' + hname + '" \n \n')
#a = """
#////////////////////////////////////////////////////////////////////////////////
#
#namespace {
#
#struct variable {
#    double v_exp0(const double& r0, const double& k,
#                 const double * p1, const double * p2 );
#                 
#    double v_exp(const double& k,
#                 const double * p1, const double * p2 );
#
#    double v_coul0(const double& r0, const double& k,
#                  const double * p1, const double * p2 );
#                  
#    double v_coul(const double& k,
#                  const double * p1, const double * p2 );
#
#    double v_gau0(const double& r0, const double& k,
#                 const double * p1, const double * p2 );
#                  
#    void grads(const double& gg, double * grd1, double * grd2,
#               const double * p1, const double * p2);
#
#    double g[3]; // diff(value, p1 - p2)
#};
#
#//----------------------------------------------------------------------------//
#
#double variable::v_gau0(const double& r0, const double& k,
#                       const double * p1, const double * p2)
#{
#    g[0] = p1[0] - p2[0];
#    g[1] = p1[1] - p2[1];
#    g[2] = p1[2] - p2[2];
#
#    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
#
#    const double exp1 = std::exp(-k*(r0 - r)*(r0 - r));
#    const double gg = 2*k*(r0 - r)*exp1/r;
#
#    g[0] *= gg;
#    g[1] *= gg;
#    g[2] *= gg;
#
#    return exp1;
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#double variable::v_exp0(const double& r0, const double& k,
#                       const double * p1, const double * p2)
#{
#    g[0] = p1[0] - p2[0];
#    g[1] = p1[1] - p2[1];
#    g[2] = p1[2] - p2[2];
#
#    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
#
#    const double exp1 = std::exp(k*(r0 - r));
#    const double gg = - k*exp1/r;
#
#    g[0] *= gg;
#    g[1] *= gg;
#    g[2] *= gg;
#
#    return exp1;
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#double variable::v_exp(const double& k,
#                       const double * p1, const double * p2)
#{
#    g[0] = p1[0] - p2[0];
#    g[1] = p1[1] - p2[1];
#    g[2] = p1[2] - p2[2];
#
#    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
#
#    const double exp1 = std::exp(k*(- r));
#    const double gg = - k*exp1/r;
#
#    g[0] *= gg;
#    g[1] *= gg;
#    g[2] *= gg;
#
#    return exp1;
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#
#double variable::v_coul(const double& k,
#                        const double * p1, const double * p2)
#{
#    g[0] = p1[0] - p2[0];
#    g[1] = p1[1] - p2[1];
#    g[2] = p1[2] - p2[2];
#
#    const double rsq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
#    const double r = std::sqrt(rsq);
#
#    const double exp1 = std::exp(k*(-r));
#    const double rinv = 1.0/r;
#    const double val = exp1*rinv;
#
#    const double gg = - (k + rinv)*val*rinv;
#
#    g[0] *= gg;
#    g[1] *= gg;
#    g[2] *= gg;
#
#    return val;
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#double variable::v_coul0(const double& r0, const double& k,
#                        const double * p1, const double * p2)
#{
#    g[0] = p1[0] - p2[0];
#    g[1] = p1[1] - p2[1];
#    g[2] = p1[2] - p2[2];
#
#    const double rsq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
#    const double r = std::sqrt(rsq);
#
#    const double exp1 = std::exp(k*(r0 - r));
#    const double rinv = 1.0/r;
#    const double val = exp1*rinv;
#
#    const double gg = - (k + rinv)*val*rinv;
#
#    g[0] *= gg;
#    g[1] *= gg;
#    g[2] *= gg;
#
#    return val;
#}
#
#//----------------------------------------------------------------------------//
#
#void variable::grads(const double& gg, double * grd1, double * grd2, 
#                     const double * p1, const double * p2) {
#    for (size_t i = 0; i < 3 ; i++) {
#        double d = gg*g[i];
#        grd1[i] += d;
#        grd2[i] -= d;
#    }
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#struct monomer {
#    double oh1[3];
#    double oh2[3];
#
#    void setup(const double* ohh,
#               const double& in_plane_g, const double& out_of_plane_g,
#               double x1[3], double x2[3]);
#
#    void grads(const double* g1, const double* g2,
#               const double& in_plane_g, const double& out_of_plane_g,
#               double* grd) const;
#};
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#void monomer::setup(const double* ohh,
#                    const double& in_plane_g, const double& out_of_plane_g,
#                    double* x1, double* x2)
#{
#    for (int i = 0; i < 3; ++i) {
#        oh1[i] = ohh[i + 3] - ohh[i];
#        oh2[i] = ohh[i + 6] - ohh[i];
#    }
#
#    const double v[3] = {
#        oh1[1]*oh2[2] - oh1[2]*oh2[1],
#        oh1[2]*oh2[0] - oh1[0]*oh2[2],
#        oh1[0]*oh2[1] - oh1[1]*oh2[0]
#    };
#
#    for (int i = 0; i < 3; ++i) {
#        const double in_plane = ohh[i] + 0.5*in_plane_g*(oh1[i] + oh2[i]);
#        const double out_of_plane = out_of_plane_g*v[i];
#
#        x1[i] = in_plane + out_of_plane;
#        x2[i] = in_plane - out_of_plane;
#    }
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#void monomer::grads(const double* g1, const double* g2,
#                    const double& in_plane_g, const double& out_of_plane_g,
#                    double* grd) const
#{
#    const double gm[3] = {
#        g1[0] - g2[0],
#        g1[1] - g2[1],
#        g1[2] - g2[2]
#    };
#
#    const double t1[3] = {
#        oh2[1]*gm[2] - oh2[2]*gm[1],
#        oh2[2]*gm[0] - oh2[0]*gm[2],
#        oh2[0]*gm[1] - oh2[1]*gm[0]
#    };
#
#    const double t2[3] = {
#        oh1[1]*gm[2] - oh1[2]*gm[1],
#        oh1[2]*gm[0] - oh1[0]*gm[2],
#        oh1[0]*gm[1] - oh1[1]*gm[0]
#    };
#
#    for (int i = 0; i < 3; ++i) {
#        const double gsum = g1[i] + g2[i];
#        const double in_plane = 0.5*in_plane_g*gsum;
#
#        const double gh1 = in_plane + out_of_plane_g*t1[i];
#        const double gh2 = in_plane - out_of_plane_g*t2[i];
#
#        grd[i + 0] += gsum - (gh1 + gh2); // O
#        grd[i + 3] += gh1; // H1
#        grd[i + 6] += gh2; // H2
#    }
#}
#
#//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
#
#//struct vsites {
#//    //void TwoParticleAverageSite() {}
#//    //void ThreeParticleAverageSite() {}
#//    void OutOfPlaneSite(const double& w12, const double& w13,
#//                        const double& wcross, const double x1[3],
#//                        const double y1[3], const double y2[3],
#//                        double vs[3]);
#//    //void LocalCoordinatesSite{}
#//};
#//
#//void vsites::OutOfPlaneSite(const double& w12,
#//                            const double& w13,
#//                            const double& wcross,
#//                            const double x1[3],
#//                            const double y1[3],
#//                            const double y2[3],
#//                            double vs[3]) {
#//    double r12[3], r13[3];
#//
#//    for (int i = 0; i < 3; ++i) {
#//        r12[i] = y1[i] - x1[i];
#//        r13[i] = y2[i] - x1[i];
#//    }
#//                            
#//    double rc[3];
#//    rc[0] = r12[1]*r13[2] - r12[2]*r13[1];
#//    rc[1] = r12[2]*r13[0] - r12[0]*r13[2];
#//    rc[2] = r12[0]*r13[1] - r12[1]*r13[0];
#//    
#//    vs[0] = x1[0] + w12 * r12[0] + w13 * r13[0] + wcross * rc[0];
#//    vs[1] = x1[1] + w12 * r12[1] + w13 * r13[1] + wcross * rc[1];
#//    vs[2] = x1[2] + w12 * r12[2] + w13 * r13[2] + wcross * rc[2];
#//}
#
#} // namespace
#
#////////////////////////////////////////////////////////////////////////////////
#
#namespace x2b_""" + mon1 + "_" + mon2 + "_deg" + str(degree) + """ {
#
#//----------------------------------------------------------------------------//
#
#x2b_""" + mon1 + "_" + mon2 + """_v1x::x2b_""" + mon1 + "_" + mon2 + """_v1x(std::string mon1, std::string mon2) {
#
#    // =====>> SECTION CONSTRUCTOR <<=====
#    // =>> PASTE RIGHT BELOW THIS LINE <==
#
#}
#
#//----------------------------------------------------------------------------//
#
#double x2b_""" + mon1 + "_" + mon2 + """_v1x::f_switch(const double& r, double& g) const
#{
#    if (r > m_r2f) {
#        g = 0.0;
#        return 0.0;
#    } else if (r > m_r2i) {
#        const double t1 = M_PI/(m_r2f - m_r2i);
#        const double x = (r - m_r2i)*t1;
#        g = - std::sin(x)*t1/2.0;
#        return (1.0 + std::cos(x))/2.0;
#    } else {
#        g = 0.0;
#        return 1.0;
#    }
#}
#
#//----------------------------------------------------------------------------//
#
#double x2b_""" + mon1 + "_" + mon2 + """_v1x::eval(const double* xyz1, const double* xyz2, const size_t ndim) const
#{
#
#    std::vector<double> energies(ndim,0.0);
#
#    for (size_t j = 0; j < ndim; j++) {
#        double mon1[""" + str(nat1*3) + """];
#        double mon2[""" + str(nat2*3) + """];
#
#        std::copy(xyz1 + j * """ + str(nat1*3) + """, xyz1 + (j+1) * """ + str(nat1*3) + """, mon1);
#        std::copy(xyz2 + j * """ + str(nat2*3) + """, xyz2 + (j+1) * """ + str(nat2*3) + """, mon2);
#
#        // Right now it assumes 1st atom of each monomer
#        const double d12[3] = {mon1[0] -  mon2[0],
#                               mon1[1] -  mon2[1],
#                               mon1[2] -  mon2[2]};
#    
#        const double r12sq = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
#        const double r12 = std::sqrt(r12sq);
#    
#        if (r12 > m_r2f)
#            continue;
#    
#        double xcrd[""" + str(3*(nat1 + nat2)) + """]; // coordinates of real sites ONLY
#    
#        std::copy(mon1, mon1 + """ + str(3*(nat1)) + """, xcrd);
#        std::copy(mon2, mon2 + """ + str(3*(nat2)) + """, xcrd + """ + str(3*(nat1)) + """);
#        
#        double v[""" + str(nvars) + """];
#        
#        double sw = 0.0;
#        double gsw = 0.0;
#    
#"""
#ff.write(a)
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_a),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_a[i+1])):
#        if not types_a[i] in vsites:
#            ff.write('        const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xcrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_a),2):
#    n = 1
#    for j in range(int(types_a[i+1])):
#        if types_a[i] in vsites:
#            ff.write('        double ' + types_a[i] + '_' + str(n) + '_a[3]' + ';\n')
#            n = n + 1
#ff.write('\n')
#
## loops over each type of atom in the input
#for i in range(0,len(types_b),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_b[i+1])):
#        if not types_b[i] in vsites:
#            ff.write('        const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xcrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_b),2):
#    n = 1
#    for j in range(int(types_b[i+1])):
#        if types_b[i] in vsites:
#            ff.write('        double ' + types_b[i] + '_' + str(n) + '_b[3]' + ';\n')
#            n = n + 1
#ff.write('\n')
#
#if use_lonepairs[0] == 0 and use_lonepairs[1] == 0:
#    a = """    
#        variable vr[""" + str(nvars) + """];
#    
#"""
#elif use_lonepairs[0] != 0:
#    a = """
#        // vsites virt;
#        double w12 =     -9.721486914088159e-02;  //from MBpol
#        double w13 =     -9.721486914088159e-02;
#        double wcross =   9.859272078406150e-02;
#
#        monomer m;
#        
#        m.setup(""" + types_a[0] + '_1_a' + """, w12, wcross,
#                 """ + types_a[4] + '_1_a' + """, """ + types_a[4] + '_2_a' + """);
#                        
#        variable vr[""" + str(nvars) + """];
#    
#"""
#elif use_lonepairs[1] != 0:
#    a = """
#        // vsites virt;
#        double w12 =     -9.721486914088159e-02;  //from MBpol
#        double w13 =     -9.721486914088159e-02;
#        double wcross =   9.859272078406150e-02;
#
#        monomer m;
#        
#        m.setup(""" + types_b[0] + '_1_b' + """, w12, wcross,
#                 """ + types_b[4] + '_1_b' + """, """ + types_b[4] + '_2_b' + """);
#                        
#        variable vr[""" + str(nvars) + """];
#    
#"""
#ff.write(a)
#
#nv = 0
## Intramolecular distances:
#for index, variable in enumerate(variables):
#    if variable[4].startswith("x-intra-"):
#        atom1 = variable[0][0]
#        atom2 = variable[2][0]
#
#        try:
#            atom1_index = int(variable[0][1:])
#        except ValueError:
#            atom1_index = 1
#
#        try:
#            atom2_index = int(variable[2][1:])
#        except ValueError:
#            atom2_index = 1
#
#        atom1_fragment = variable[1]
#        atom2_fragment = variable[3]
#
#        atom1_name = "{}_{}_{}".format(atom1, atom1_index, atom1_fragment)
#        atom2_name = "{}_{}_{}".format(atom2, atom2_index, atom2_fragment)
#
#        sorted_atoms = "".join(sorted([atom1, atom2]))
#
#        if var_intra == 'exp0' or var_intra == 'coul0' or var_intra == 'gau0':
#            arguments = '(m_d_intra_' + sorted_atoms + ', m_k_intra_' + sorted_atoms
#        else:
#            arguments = '(m_k_intra_' + sorted_atoms
#
#        ff.write('        v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var_intra + arguments + ', ' + atom1_name + ', ' + atom2_name + ');\n')
#
#        nv += 1
#            
#    else:
#        # inter-molecular variables
#        atom1 = variable[0][0]
#        atom2 = variable[2][0]
#
#        try:
#            atom1_index = int(variable[0][1:])
#        except ValueError:
#            atom1_index = 1
#
#        try:
#            atom2_index = int(variable[2][1:])
#        except ValueError:
#            atom2_index = 1
#
#        atom1_fragment = variable[1]
#        atom2_fragment = variable[3]
#
#        atom1_name = "{}_{}_{}".format(atom1, atom1_index, atom1_fragment)
#        atom2_name = "{}_{}_{}".format(atom2, atom2_index, atom2_fragment)
#
#        if atom1 not in vsites and atom2 not in vsites:
#            var = var_inter
#        else:
#            var = var_lp
#
#        sorted_atoms = "".join(sorted([atom1, atom2]))
#
#        if var == 'exp0' or var == 'coul0' or var == 'gau0':
#            arguments = '(m_d_' + sorted_atoms + ', m_k_' + sorted_atoms
#        else:
#            arguments = '(m_k_' + sorted_atoms
#
#        ff.write('        v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var + arguments + ', ' + atom1_name + ', ' + atom2_name + ');\n')
#
#        nv += 1
# 
#a = """     
#    
#        sw = f_switch(r12, gsw);
#        
#        energies[j] = sw*polynomial::eval(coefficients.data(), v);
#    }
#
#    double energy = 0.0;
#    for (size_t i = 0; i < ndim; i++) {
#      energy += energies[i];
#    }
#
#    return energy;
#    
#}
#
#double x2b_""" + mon1 + "_" + mon2 + """_v1x::eval(const double* xyz1, const double* xyz2, 
#                double * grad1, double * grad2, const size_t ndim) const
#{
#
#    std::vector<double> energies(ndim,0.0);
#
#    for (size_t j = 0; j < ndim; j++) {
#        double mon1[""" + str(nat1*3) + """];
#        double mon2[""" + str(nat2*3) + """];
#
#        std::copy(xyz1 + j * """ + str(nat1*3) + """, xyz1 + (j+1) * """ + str(nat1*3) + """, mon1);
#        std::copy(xyz2 + j * """ + str(nat2*3) + """, xyz2 + (j+1) * """ + str(nat2*3) + """, mon2);
#
#        // Right now it assumes 1st atom of each monomer
#        const double d12[3] = {mon1[0] -  mon2[0],
#                               mon1[1] -  mon2[1],
#                               mon1[2] -  mon2[2]};
#    
#        const double r12sq = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
#        const double r12 = std::sqrt(r12sq);
#    
#        if (r12 > m_r2f)
#            continue;
#    
#        double xcrd[""" + str(3*(nat1 + nat2)) + """]; // coordinates of real sites ONLY
#    
#        std::copy(mon1, mon1 + """ + str(3*(nat1)) + """, xcrd);
#        std::copy(mon2, mon2 + """ + str(3*(nat2)) + """, xcrd + """ + str(3*(nat1)) + """);
#        
#        double v[""" + str(nvars) + """];
#        
#        double sw = 0.0;
#        double gsw = 0.0;
#    
#"""
#ff.write(a)
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_a),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_a[i+1])):
#        if not types_a[i] in vsites:
#            ff.write('        const double* ' + types_a[i] + '_' + str(n) + '_a' + '= xcrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_a),2):
#    n = 1
#    for j in range(int(types_a[i+1])):
#        if types_a[i] in vsites:
#            ff.write('        double ' + types_a[i] + '_' + str(n) + '_a[3]' + ';\n')
#            n = n + 1
#ff.write('\n')
#
## loops over each type of atom in the input
#for i in range(0,len(types_b),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_b[i+1])):
#        if not types_b[i] in vsites:
#            ff.write('        const double* ' + types_b[i] + '_' + str(n) + '_b' + '= xcrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_b),2):
#    n = 1
#    for j in range(int(types_b[i+1])):
#        if types_b[i] in vsites:
#            ff.write('        double ' + types_b[i] + '_' + str(n) + '_b[3]' + ';\n')
#            n = n + 1
#ff.write('\n')
#
#if use_lonepairs[0] == 0 and use_lonepairs[1] == 0:
#    a = """    
#        variable vr[""" + str(nvars) + """];
#    
#"""
#elif use_lonepairs[0] != 0:
#    a = """
#        //vsites virt;
#        double w12 =     -9.721486914088159e-02;  //from MBpol
#        double w13 =     -9.721486914088159e-02;
#        double wcross =   9.859272078406150e-02;
#
#        monomer m;
#        
#        m.setup(""" + types_a[0] + '_1_a' + """, w12, wcross, 
#                 """ + types_a[4] + '_1_a' + """, """ + types_a[4] + '_2_a' + """);
#                        
#        variable vr[""" + str(nvars) + """];
#    
#"""
#elif use_lonepairs[1] != 0:
#    a = """
#        //vsites virt;
#        double w12 =     -9.721486914088159e-02;  //from MBpol
#        double w13 =     -9.721486914088159e-02;
#        double wcross =   9.859272078406150e-02;
#
#        monomer m;
#        
#        m.setup(""" + types_b[0] + '_1_b' + """, w12, wcross, 
#                 """ + types_b[4] + '_1_b' + """, """ + types_b[4] + '_2_b' + """);
#                        
#        variable vr[""" + str(nvars) + """];
#    
#"""
#ff.write(a)
#
#nv = 0
## Intramolecular distances:
#for index, variable in enumerate(variables):
#    if variable[4].startswith("x-intra-"):
#        atom1 = variable[0][0]
#        atom2 = variable[2][0]
#
#        try:
#            atom1_index = int(variable[0][1:])
#        except ValueError:
#            atom1_index = 1
#
#        try:
#            atom2_index = int(variable[2][1:])
#        except ValueError:
#            atom2_index = 1
#
#        atom1_fragment = variable[1]
#        atom2_fragment = variable[3]
#
#        atom1_name = "{}_{}_{}".format(atom1, atom1_index, atom1_fragment)
#        atom2_name = "{}_{}_{}".format(atom2, atom2_index, atom2_fragment)
#
#        sorted_atoms = "".join(sorted([atom1, atom2]))
#
#        if var_intra == 'exp0' or var_intra == 'coul0' or var_intra == 'gau0':
#            arguments = '(m_d_intra_' + sorted_atoms + ', m_k_intra_' + sorted_atoms
#        else:
#            arguments = '(m_k_intra_' + sorted_atoms
#
#        ff.write('        v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var_intra + arguments + ', ' + atom1_name + ', ' + atom2_name + ');\n')
#
#        nv += 1
#            
#    else:
#        # inter-molecular variables
#        atom1 = variable[0][0]
#        atom2 = variable[2][0]
#
#        try:
#            atom1_index = int(variable[0][1:])
#        except ValueError:
#            atom1_index = 1
#
#        try:
#            atom2_index = int(variable[2][1:])
#        except ValueError:
#            atom2_index = 1
#
#        atom1_fragment = variable[1]
#        atom2_fragment = variable[3]
#
#        atom1_name = "{}_{}_{}".format(atom1, atom1_index, atom1_fragment)
#        atom2_name = "{}_{}_{}".format(atom2, atom2_index, atom2_fragment)
#
#        if atom1 not in vsites and atom2 not in vsites:
#            var = var_inter
#        else:
#            var = var_lp
#
#        sorted_atoms = "".join(sorted([atom1, atom2]))
#
#        if var == 'exp0' or var == 'coul0' or var == 'gau0':
#            arguments = '(m_d_' + sorted_atoms + ', m_k_' + sorted_atoms
#        else:
#            arguments = '(m_k_' + sorted_atoms
#
#        ff.write('        v[' + str(nv) + ']  = vr[' + str(nv) + '].v_' + var + arguments + ', ' + atom1_name + ', ' + atom2_name + ');\n')
#
#        nv += 1
#
#a = """     
#    
#        double g[""" + str(nvars) + """];
#
#        // the switch
#        sw = f_switch(r12, gsw);
#        
#        energies[j] = sw*polynomial::eval(coefficients.data(), v, g);
#        
#        double xgrd[""" + str(3*(len(atom_list_a) + len(atom_list_b))) + """];
#        std::fill(xgrd, xgrd + """ + str(3*(len(atom_list_a) + len(atom_list_b))) + """, 0.0);
#
#""" 
#ff.write(a)   
#
#nc = 0
## loops over each type of atom in the input
#for i in range(0,len(types_a),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_a[i+1])):
#        if not types_a[i] in vsites:
#            ff.write('        double* ' + types_a[i] + '_' + str(n) + '_a_g' + '= xgrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
## loops over each type of atom in the input
#for i in range(0,len(types_b),2):
#    n = 1
#    # loops over each atom of that type
#    for j in range(int(types_b[i+1])):
#        if not types_b[i] in vsites:
#            ff.write('        double* ' + types_b[i] + '_' + str(n) + '_b_g' + '= xgrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#            nc = nc + 1
#ff.write('\n')
#
#for i in range(0,len(types_a),2):
#    n = 1
#    for j in range(int(types_a[i+1])):
#        if types_a[i] in vsites:
#            ff.write('        double* ' + types_a[i] + '_' + str(n) + '_a_g' + '= xgrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#ff.write('\n')
#
#for i in range(0,len(types_b),2):
#    n = 1
#    for j in range(int(types_b[i+1])):
#        if types_b[i] in vsites:
#            ff.write('        double* ' + types_b[i] + '_' + str(n) + '_b_g' + '= xgrd + ' + str(3 * nc) + ';\n')
#            n = n + 1
#ff.write('\n')
#
#nv = 0
## Intramolecular distances:
#for index, variable in enumerate(variables):
#    atom1 = variable[0][0]
#    atom2 = variable[2][0]
#
#    try:
#        atom1_index = int(variable[0][1:])
#    except ValueError:
#        atom1_index = 1
#
#    try:
#        atom2_index = int(variable[2][1:])
#    except ValueError:
#        atom2_index = 1
#
#    atom1_fragment = variable[1]
#    atom2_fragment = variable[3]
#
#    atom1_name = "{}_{}_{}".format(atom1, atom1_index, atom1_fragment)
#    atom2_name = "{}_{}_{}".format(atom2, atom2_index, atom2_fragment)
#
#    ff.write('        vr[' + str(nv) + '].grads(g[' + str(nv) + '], ' + atom1_name + '_g, ' + atom2_name + '_g, ' + atom1_name + ', ' + atom2_name + ');\n')
#
#    nv += 1
#            
#a = """
#
#    // ##DEFINE HERE## the redistribution of the gradients
#    
#"""
#ff.write(a)
#a = ""
#if use_lonepairs[0] != 0:
#    a = """
#        m.grads(""" + types_a[4] + '_1_a_g, ' + types_a[4] + """_2_a_g, 
#                 w12, wcross, """ + types_a[0] + """_1_a_g);
#    """
#elif use_lonepairs[1] != 0:
#    a = """
#        m.grads(""" + types_b[4] + '_1_b_g, ' + types_b[4] + """_2_b_g, 
#                 w12, wcross, """ + types_b[0] + """_1_b_g);
#"""
#ff.write(a)    
#
#a = """
#    
#        for (int i = 0; i < """ + str(3*nat1) + """; ++i) {
#            grad1[i + j*""" + str(3*nat1) + """] += sw*xgrd[i];
#        }
#
#        for (int i = 0; i < """ + str(3*nat2) + """; ++i) {
#            grad2[i + j*""" + str(3*nat2) + """] += sw*xgrd[i + """ + str(3*nat1) + """];
#        }
#
#        // gradient of the switch
#
#        gsw *= energies[j]/r12;
#        for (int i = 0; i < 3; ++i) {
#            const double d = gsw*d12[i];
#            grad1[i + j*""" + str(3*nat1) + """] += d;
#            grad2[i + j*""" + str(3*nat2) + """] -= d;
#        }
#
#    }
#
#    double energy = 0.0;
#    for (size_t i = 0; i < ndim; i++) {
#      energy += energies[i];
#    }
#
#    return energy;
#}
#
#} // namespace x2b_""" + mon1 + "_" + mon2 + "_deg" + str(degree) + """
#
#////////////////////////////////////////////////////////////////////////////////
#"""
#ff.write(a)
#ff.close()
#