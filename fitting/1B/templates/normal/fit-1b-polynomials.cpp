#include <cmath>
#include <cassert>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include <gsl/gsl_multimin.h>

#define RIDGE_REGRESSION 1

#ifdef RIDGE_REGRESSION
#   include "rwlsq.h"
#else
#   include "wlsq.h"
#endif

#include "x1b-v1.h"
#include "training-set.h"
#include "fit-utils.h"

#define dont_be_VERBOSE yes

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

typedef h3o::x1b_v1 model_type;
static model_type model;

#ifdef RIDGE_REGRESSION
const double alpha = 0.0001;
#endif

const double E_range = 50.0; // kcal/mol

//----------------------------------------------------------------------------//

static std::vector<h3o::monomer> training_set;
static double*                  ts_weights = 0;

////////////////////////////////////////////////////////////////////////////////

namespace linear {

//----------------------------------------------------------------------------//

static double* A(0);
static double* y(0);
static double* params(0);

//----------------------------------------------------------------------------//

void allocate()
{
    A = new double[training_set.size()*model_type::nparams()];
    y = new double[training_set.size()];

    params = new double[model_type::nparams()];
}

//----------------------------------------------------------------------------//

double compute_chisq(const gsl_vector* X, void* unused)
{
    model.set_nonlinear_parameters(X->data);

    if (model.nonlinear_parameters_out_of_range())
        return 1.0e+6;

    for (size_t n = 0; n < training_set.size(); ++n) {
        double mmm[model_type::nparams()];
        y[n] = training_set[n].energy_onebody
             - model.basis(training_set[n].xyz, mmm);
        for (size_t p = 0; p < model_type::nparams(); ++p)
            A[p + n*model_type::nparams()] = mmm[p];
    }

#   ifdef VERBOSE
    std::cout << "=== calling wlsq::solve() ["
              << kit::wlsq::implementation()
              << "] ===" << std::endl;
#   endif

    double chisq;

#   ifdef RIDGE_REGRESSION
    double penaltysq;
//    std::cout << "inside ridge regression. calling solve.\n";
    kit::rwlsq::solve(training_set.size(), model_type::nparams(),
                      A, y, ts_weights, alpha, params, chisq, penaltysq);

    std::cout << "<#> chisq = " << chisq
              << " : penaltysq = " << penaltysq
              << std::endl;
#   else
    int rank;
    kit::wlsq::solve(training_set.size(), model_type::nparams(),
                     A, y, ts_weights, params, chisq, rank);
    std::cout << "<#> chisq = " << chisq
              << " : rank = " << rank
              << std::endl;
#   endif

#   ifdef VERBOSE
    std::cout << "\n--> chisq = " << chisq
              << "\n-->  rank = " << rank
              << '\n' << std::endl;
#   endif

//    delete[] A;
//    delete[] y;
//    delete[] params;
    return chisq;
}

//----------------------------------------------------------------------------//

} // namespace linear

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    //
    // initialize the model and load training sets
    //

    if (argc < 2) {
        std::cerr << "usage: fit-x1b-v1 ts1 ..."
                  << std::endl;
        return 0;
    }

    ++argv; --argc;

    std::cout << "\n<><><> model type = '" << model.name() << "'\n";
#   ifdef RIDGE_REGRESSION
    std::cout << "<> using ridge regression with alpha = "  << alpha << std::endl;
#   endif

    std::cout << "\n<><><> model type = '" << model.name() << "'\n";

    const double x0[] = {
        1.0,
        1.0,
	2.0,
        1.0,
	2.0
    };
    model.set_nonlinear_parameters(x0);
//    model.set_r3i_r3f(0.0, 4.5);
//    std::cout << "non-linear parameters set \n";

    const char fn[] = "fit-h3o-1b-v1-initial.dat";
    std::ofstream ofs(fn);
    std::cout << "\n>> dumping initial model as '" << fn << "' >>\n\n";
    model.as_poly_dat(ofs);
    std::cout << "initial value written\n" ;
    ofs.close();

    try {
        while (argc-- != 0) {
            size_t nt = h3o::load_monomers(*argv, training_set);
            std::cout << "'" << *(argv++) << "' : "
                      << nt << " monomers" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    //
    // assign weights
    //

    ts_weights = new double[training_set.size()];
//    std::fill(ts_weights, ts_weights + training_set.size(), 1.0);

    double E_min(0.0), E_max(0.0);
    for (size_t n = 0; n < training_set.size(); ++n) {
        if (E_min > training_set[n].energy_onebody)
            E_min = training_set[n].energy_onebody;
        if (E_max < training_set[n].energy_onebody)
            E_max = training_set[n].energy_onebody;
   }

    size_t N_eff;

    h3o::setup_weights(training_set, E_range, E_min,
                       ts_weights, N_eff);

    std::cout << "\n>>   E_min = " << E_min << " kcal/mol"
                 "\n>> E_range = " << E_range << " kcal/mol"
                 "\n\n>> training set size = " << training_set.size()
              << "\n>>    effective size = " << N_eff << '\n'
              << std::endl;


//    std::cout << ">> E_min/E_max = " << E_min << " / " << E_max << std::endl;

//    std::cout << ">> " << training_set.size() << " monomers \n";

    //
    // set things up
    //

    linear::allocate();

    gsl_vector* x = gsl_vector_alloc(model_type::num_nonlinear_params);
    model.get_nonlinear_parameters(x->data);

    gsl_multimin_function chisq_func;

    chisq_func.n = model_type::num_nonlinear_params;
    chisq_func.f = linear::compute_chisq;

    gsl_vector* ss = gsl_vector_alloc(model_type::num_nonlinear_params);
    gsl_vector_set_all(ss, 0.3);

    std::cout << "\n<> initial simplex sides:\n";
    for (size_t n = 0; n < model_type::num_nonlinear_params; ++n)
        std::cout << n << " : " << ss->data[n] << "\n";

    gsl_multimin_fminimizer* s =
        gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2rand,
                                      model_type::num_nonlinear_params);

    gsl_multimin_fminimizer_set(s, &chisq_func, x, ss);

    int iter(0), status(0);

    do {
        ++iter;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status)
             break;

        const double size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-4);  //TODO convergence limit

        if (status == GSL_SUCCESS)
            std::cout << "!!! converged !!!" << std::endl;

        std::cout << iter << ' ' << s->fval << ' ' << size << std::endl;

        if (iter > 0 && iter%10 == 0) {
            std::cout << "\n<> solution:\n";
            for (size_t n = 0; n <  model_type::num_nonlinear_params; ++n)
                std::cout << n << " : " << s->x->data[n] << "\n";
            std::cout << "<>" << std::endl;
        }
    } while (status == GSL_CONTINUE && iter < 2500);

    model.set_nonlinear_parameters(s->x->data);
    linear::compute_chisq(s->x, 0);
    model.set(linear::params);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    //
    // report
    //


    std::ofstream correlation_file;
    correlation_file.open ("correlation.dat");
    double err_L2(0), err_wL2(0), err_Linf(0);
    double err_L2_lo(0), err_Linf_lo(0), nlo(0);

    for (size_t i = 0; i < training_set.size(); ++i) {

        const double E_model = model.value(training_set[i].xyz);
        const double delta = E_model - training_set[i].energy_onebody;
        if (std::abs(delta) > err_Linf)
            err_Linf = std::abs(delta);

        correlation_file // <<  i+1   << "   "
          <<  training_set[i].energy_onebody  << "    "
          << E_model  << "    \n" ;
        //  <<  delta*delta  << "    \n" ;

        err_L2 += delta*delta;
        err_wL2 += ts_weights[i]*delta*delta;

//        if (training_set[i].energy_total - E_min < E_range) {
//            nlo += 1.0;
//            err_L2_lo += delta*delta;
//            if (std::abs(delta) > err_Linf_lo)
//                err_Linf_lo = std::abs(delta);
//        }
    }
    
    correlation_file.close();

    // RMSD
    err_L2 /= training_set.size();

    // weighted-RMSD 
    err_wL2 /= training_set.size();

    // RMSD for low-energy configurations
    //      low energy given by E_range
    err_L2_lo /= nlo;

    std::cout << "      err[L2] = " << std::sqrt(err_L2) << "\n"
              << "     err[wL2] = " << std::sqrt(err_wL2) << "\n"
              << "    err[Linf] = " << err_Linf << "\n"
              << std::endl;

    //
    //  save
    //

    {
        const char fn_fin[] = "fit-x1b-v1.dat";
        std::ofstream ofs(fn_fin);
        std::cout << "\n>> saving as '" << fn_fin << "'\n";
        model.as_poly_dat(ofs);
    }

//    delete[] ts_weights;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
