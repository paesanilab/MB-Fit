
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <chrono>

#include <gsl/gsl_multimin.h>

#define RIDGE_REGRESSION 1

#ifdef RIDGE_REGRESSION
#   include "rwlsq.h"
#else
#   include "wlsq.h"
#endif

#include "mon1.h"
#include "mon2.h"
#include "fit-utils.h"
#include "training_set.h"
#include "x2b_A1B2Z2_D1E2_v1.h"
#include "electrostatics.h"
#include "coulomb.h"
#include "dispersion.h"


#define dont_be_VERBOSE yes

//#define DEBUG 

#ifdef DEBUG
#define PR(x) std::cout << #x << ": " << (x) << std::endl;
#else
#define PR(x)
#endif /* DEBUG */

using namespace std;

namespace {

static x2b_A1B2Z2_D1E2::x2b_A1B2Z2_D1E2_v1x model;


#ifdef RIDGE_REGRESSION
const double alpha = 0.0005;
#endif

// ##DEFINE HERE## energy range
const double E_range = 30.0; // kcal/mol

//----------------------------------------------------------------------------//

static std::vector<tset::dimer> training_set;
static std::vector<double> elec_e;
static std::vector<double> disp_e;
static double* ts_weights = 0;

namespace linear {

//----------------------------------------------------------------------------//

static double* A(0);
static double* y(0);
static double* params(0);

//----------------------------------------------------------------------------//

void allocate()
{
    A = new double[training_set.size()*model.nparams()];
    y = new double[training_set.size()];

    params = new double[model.nparams()];
}

//----------------------------------------------------------------------------//

double compute_chisq(const gsl_vector* X, void* unused)
{
    model.set_nonlinear_parameters(X->data);

    if (model.nonlinear_parameters_out_of_range())
        return 1.0e+6;

    PR(model_type::nparams())
    for (size_t n = 0; n < training_set.size(); ++n) {
        double mmm[model.nparams()];
        y[n] = training_set[n].energy_twobody
               - model.basis(training_set[n].xyz, mmm);
        for (size_t p = 0; p < model.nparams(); ++p) {
            A[p + n*model.nparams()] = mmm[p];
        }
    }

#   ifdef VERBOSE
    std::cout << "=== calling wlsq::solve() ["
              << kit::wlsq::implementation()
              << "] ===" << std::endl;
#   endif

    double chisq;

//    try {
#     ifdef RIDGE_REGRESSION
      double penaltysq;
      kit::rwlsq::solve(training_set.size(), model.nparams(),
                        A, y, ts_weights, alpha, params, chisq, penaltysq);

      std::cout << "<#> chisq = " << chisq
                << " : penaltysq = " << penaltysq
                << std::endl;
#     else
      int rank;
      kit::wlsq::solve(training_set.size(), model_type::nparams(),
                       A, y, ts_weights, params, chisq, rank);
      std::cout << "<#> chisq = " << chisq
                << " : rank = " << rank
                << std::endl;
#     endif

      if (!gsl_finite (chisq)) {
        return 10000000.0;
      }
//    } catch (const std::exception& e) {
//      return 1000000.0;
//    }
#   ifdef VERBOSE
    std::cout << "\n--> chisq = " << chisq
              << "\n-->  rank = " << rank
              << '\n' << std::endl;
#   endif

    return chisq;
}

//----------------------------------------------------------------------------//

} // namespace linear

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: fit-2b ts1 ..."
                  << std::endl;
        return 0;
    }
      
    long long int duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::system_clock::now().time_since_epoch()).count();

    srand(duration);

    double x0[20];
      x0[0] = ((double) rand() / (RAND_MAX)) * 3.0 + 0.0;
      x0[1] = ((double) rand() / (RAND_MAX)) * 2.0 + 0.0;
      x0[2] = ((double) rand() / (RAND_MAX)) * 3.0 + 0.0;
      x0[3] = ((double) rand() / (RAND_MAX)) * 2.0 + 0.0;
      x0[4] = ((double) rand() / (RAND_MAX)) * 3.0 + 0.0;
      x0[5] = ((double) rand() / (RAND_MAX)) * 2.0 + 0.0;
      x0[6] = ((double) rand() / (RAND_MAX)) * 3.0 + 0.0;
      x0[7] = ((double) rand() / (RAND_MAX)) * 2.0 + 0.0;
      x0[8] = ((double) rand() / (RAND_MAX)) * 7.0 + 0.0;
      x0[9] = ((double) rand() / (RAND_MAX)) * 3.0 + 0.0;
      x0[10] = ((double) rand() / (RAND_MAX)) * 7.0 + 0.0;
      x0[11] = ((double) rand() / (RAND_MAX)) * 3.0 + 0.0;
      x0[12] = ((double) rand() / (RAND_MAX)) * 7.0 + 0.0;
      x0[13] = ((double) rand() / (RAND_MAX)) * 3.0 + 0.0;
      x0[14] = ((double) rand() / (RAND_MAX)) * 7.0 + 0.0;
      x0[15] = ((double) rand() / (RAND_MAX)) * 3.0 + 0.0;
      x0[16] = ((double) rand() / (RAND_MAX)) * 7.0 + 0.0;
      x0[17] = ((double) rand() / (RAND_MAX)) * 3.0 + 0.0;
      x0[18] = ((double) rand() / (RAND_MAX)) * 7.0 + 0.0;
      x0[19] = ((double) rand() / (RAND_MAX)) * 3.0 + 0.0;

      

    model.set_nonlinear_parameters(x0);

    #   ifdef RIDGE_REGRESSION
    std::cout << "<> using ridge regression with alpha = "
              << alpha << std::endl;
#   endif

    std::cout << "\n<><><> model type = '" << model.name() << "'\n";

    {
        const char fn[] = "fit-2b-initial.cdl";
        std::ofstream ofs(fn);
        std::cout << "\n>> dumping initial model as '" << fn << "' >>\n\n";
        model.as_cdl(ofs);
    }

    ++argv;

    try {
        while (--argc != 0) {
            size_t nd = tset::load_dimers(*argv, training_set);
            std::cout << "'" << *(argv++) << "' : "
                      << nd << " dimers" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    //
    // assign weights
    //

    ts_weights = new double[training_set.size()];

    double E_min;
    size_t N_eff;

    tset::setup_weights(training_set, E_range, E_min,
                       ts_weights, N_eff, true);

    std::cout << "\n>>   E_min = " << E_min << " kcal/mol"
                 "\n>> E_range = " << E_range << " kcal/mol"
                 "\n\n>> training set size = " << training_set.size()
              << "\n>>    effective size = " << N_eff << '\n'
              << std::endl;


    // electrostatics
    double tb_ref[training_set.size()];
    for (size_t n = 0; n < training_set.size(); n++) {
      // Saving reference 2b energies  
      tb_ref[n] = training_set[n].energy_twobody ;

      x::mon1 m1(training_set[n].xyz);
      x::mon2 m2(training_set[n].xyz + 3*m1.get_realsites());

      int system_nsites = m1.get_nsites() + m2.get_nsites();
      int * system_is_w;
      double* system_sitecrds;
      double* system_charge;
      double* system_polfac;
      double* system_pol;

      system_sitecrds = new double [system_nsites*3];
      system_charge = new double [system_nsites];
      system_polfac = new double [system_nsites];
      system_pol = new double [system_nsites]; //allocates memory to the pointers
      system_is_w = new int[system_nsites];
      
      std::fill(system_is_w, system_is_w + m1.get_nsites(), m1.is_w);
      std::fill(system_is_w + m1.get_nsites(), system_is_w + system_nsites, m2.is_w);

      std::copy(m1.get_sitecrds(), m1.get_sitecrds() + 3 * m1.get_nsites(),
                system_sitecrds);
      std::copy(m2.get_sitecrds(), m2.get_sitecrds() + 3 * m2.get_nsites(),
                system_sitecrds + 3 * m1.get_nsites());

      std::copy(m1.get_charges(), m1.get_charges() + m1.get_nsites(),
                system_charge);
      std::copy(m2.get_charges(), m2.get_charges() + m2.get_nsites(),
                system_charge + m1.get_nsites());

      std::copy(m1.get_pol(), m1.get_pol() + m1.get_nsites(),
                system_pol);
      std::copy(m2.get_pol(), m2.get_pol() + m2.get_nsites(),
                system_pol + m1.get_nsites());

      std::copy(m1.get_polfacs(), m1.get_polfacs() + m1.get_nsites(),
                system_polfac);
      std::copy(m2.get_polfacs(), m2.get_polfacs() + m2.get_nsites(),
                system_polfac + m1.get_nsites());


      excluded_set_type exclude12;
      excluded_set_type exclude13;
      excluded_set_type exclude14;

      for (auto i = m1.get_begin_12(); i != m1.get_end_12(); i++) {
        exclude12.insert(*i);
      }
      for (auto i = m2.get_begin_12(); i != m2.get_end_12(); i++) {
        std::pair<size_t,size_t> p =
                      std::make_pair(i->first + m1.get_nsites() ,
                                     i->second + m1.get_nsites());
        exclude12.insert(p);
      }

      for (auto i = m1.get_begin_13(); i != m1.get_end_13(); i++) {
        exclude13.insert(*i);
      }
      for (auto i = m2.get_begin_13(); i != m2.get_end_13(); i++) {
        std::pair<size_t,size_t> p =
                      std::make_pair(i->first + m1.get_nsites() ,
                                     i->second + m1.get_nsites());
        exclude13.insert(p);
      }

      for (auto i = m1.get_begin_14(); i != m1.get_end_14(); i++) {
        exclude14.insert(*i);
      }
      for (auto i = m2.get_begin_14(); i != m2.get_end_14(); i++) {
        std::pair<size_t,size_t> p =
                      std::make_pair(i->first + m1.get_nsites() ,
                                     i->second + m1.get_nsites());
        exclude14.insert(p);
      }

      ttm::electrostatics m_electrostatics;

      ttm::smear_ttm4x smr; 
      smr.m_aDD_intra_12 = 0.3;
      smr.m_aDD_intra_13 = 0.3;
      smr.m_aDD_intra_14 = 0.055;

      double ener = m_electrostatics(system_nsites, system_charge, system_polfac, system_pol,
                                    system_sitecrds, exclude12, exclude13, exclude14, 
                                    system_is_w, smr, 0);

      // Take out electrostatic energy:
      elec_e.push_back(ener);
      training_set[n].energy_twobody -= ener;
      // std::cerr << "Conf " << n << " : Elec= " << ener ;

      // Now need to take out dispersion
      x2b_disp disp(m1.get_sitecrds(), m2.get_sitecrds(), m1.get_realsites(), m2.get_realsites());
      ener = disp.get_dispersion();
      disp_e.push_back(ener);
      training_set[n].energy_twobody -= ener ;

      // std::cerr << " , Disp= " << ener << std::endl;
      delete[] system_sitecrds ;
      delete[] system_charge  ;
      delete[] system_polfac  ;
      delete[] system_pol ;
      delete[] system_is_w;

    }

        linear::allocate();

    gsl_vector* x = gsl_vector_alloc(model.get_num_nonlinear_params());
    model.get_nonlinear_parameters(x->data);

    gsl_multimin_function chisq_func;

    chisq_func.n = model.get_num_nonlinear_params();
    chisq_func.f = linear::compute_chisq;

    gsl_vector* ss = gsl_vector_alloc(model.get_num_nonlinear_params());
    gsl_vector_set_all(ss, 0.1);

    std::cout << "\n<> initial simplex sides:\n";
    for (size_t n = 0; n < model.get_num_nonlinear_params(); ++n)
        std::cout << n << " : " << ss->data[n] << "\n";

    gsl_multimin_fminimizer* s =
        gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2rand,
                                      model.get_num_nonlinear_params());

    gsl_multimin_fminimizer_set(s, &chisq_func, x, ss);

    int iter(0), status(0);
    //
    // Main optimization loop
    //

    do {
        ++iter;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
             break;

        const double size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-4);
        //status = gsl_multimin_test_size(size, 1e-3); //changed it to test convergence 

        if (status == GSL_SUCCESS)
            std::cout << "!!! converged !!!" << std::endl;

        std::cout << iter << ' ' << s->fval << ' ' << size << std::endl;
        if (iter > 0 && iter%10 == 0) {
            std::cout << "\n<> solution:\n";
            for (size_t n = 0; n <  model.get_num_nonlinear_params(); ++n)
                std::cout << n << " : " << s->x->data[n] << "\n";
            std::cout << "<>" << std::endl;
        }
    } while (status == GSL_CONTINUE && iter < 5000);

    model.set_nonlinear_parameters(s->x->data);
    linear::compute_chisq(s->x, 0);
    model.set(linear::params);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    //
    // report
    //

    for (size_t i = 0; i < training_set.size(); ++i)
      training_set[i].energy_twobody += elec_e[i] + disp_e[i];

    ofstream correlation_file;
    correlation_file.open ("correlation.dat");
    double err_L2(0), err_wL2(0), err_Linf(0);
    double err_L2_lo(0), err_Linf_lo(0), nlo(0);

    for (size_t i = 0; i < training_set.size(); ++i) {

        const double E_model = model(training_set[i].xyz) + elec_e[i] + disp_e[i];
        const double delta = E_model - tb_ref[i];
        if (std::abs(delta) > err_Linf)
            err_Linf = std::abs(delta);

  correlation_file <<  i+1   << "   "
    << E_model  << "   "
    <<  tb_ref[i]  << "    "
          <<  delta*delta  << "    \n" ;

        err_L2 += delta*delta;
        err_wL2 += ts_weights[i]*delta*delta;

        if (training_set[i].energy_total - E_min < E_range) {
            nlo += 1.0;
            err_L2_lo += delta*delta;
            if (std::abs(delta) > err_Linf_lo)
                err_Linf_lo = std::abs(delta);
        }
    }

    correlation_file.close();

    err_L2 /= training_set.size();
    err_wL2 /= training_set.size();

    err_L2_lo /= nlo;

    std::cout << "      err[L2] = " << std::sqrt(err_L2) << "    #rmsd of full ts\n"
              << "     err[wL2] = " << std::sqrt(err_wL2) << "   #weighted rmsd of full ts\n"
              << "    err[Linf] = " << err_Linf << "   #highest error in full ts\n"
              << "  err[L2,low] = " << std::sqrt(err_L2_lo) << "   #rmsd of low-energy ts \n"
              << "err[Linf,low] = " << err_Linf_lo << "   #highest error in low-energy ts "
              << std::endl;

    //
    //  save
    //

    {
        const char fn[] = "fit-2b.cdl";
        std::ofstream ofs(fn);
        std::cout << "\n>> saving as '" << fn << "'\n";
        model.as_cdl(ofs);
    }


    return 0;
}

