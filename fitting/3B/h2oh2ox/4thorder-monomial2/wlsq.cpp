#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>

#define USE_MKL

#ifdef USE_MKL
#  include <mkl_lapack.h>
#  undef USE_GSL
#else
#  define USE_GSL 1
#  include <gsl/gsl_blas.h>
#  include <gsl/gsl_vector.h>
#  include <gsl/gsl_multifit.h>
#endif // USE_MKL

#include "wlsq.h"

//mkl_set_num_threads (2);
////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

#ifdef USE_MKL
static const char implementation_string[] = "MKL";
#endif

////////////////////////////////////////////////////////////////////////////////

#ifdef USE_GSL
static const char implementation_string[] = "GSL";
#endif

////////////////////////////////////////////////////////////////////////////////

#ifdef USE_MKL
void query_workspace_sizes(int M, int N, int& lwork, int& liwork)
{
    int one(1), info, rank;
    double work;

// dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)

    lwork = -1;
    dgelsd(&M, &N, &one,
           0, &M, 0, &M, 0, 0, &rank,
           &work, &lwork, &liwork, &info);

    assert(info == 0);
    lwork = int(work);
}
#endif // USE_MKL

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace kit {

////////////////////////////////////////////////////////////////////////////////

const char* wlsq::implementation()
{
    return implementation_string;
}

////////////////////////////////////////////////////////////////////////////////

#ifdef USE_MKL
void wlsq::solve(int n_samples, int n_parameters,
     const double* A, // n_samples x n_parameters (row-major)
     const double* Y, // n_samples
     const double* w, // n_samples
           double* x, // n_parameters,
           double& chisq, int& rank)
{
    assert(n_samples > n_parameters);

    int lwork, liwork;
    query_workspace_sizes(n_samples, n_parameters, lwork, liwork);

    const size_t ndw = n_samples*n_parameters
                     + n_samples + lwork + n_parameters
                     + liwork*sizeof(double)/sizeof(int);

    double* Acopy = new double[ndw];
    double* b = Acopy + n_samples*n_parameters;
    double* work = b + n_samples;
    double* S = work + lwork;
    int* iwork = reinterpret_cast<int*>(S + n_parameters);

    // scale A according to the weights && populate b

    for (int i = 0; i < n_samples; ++i) {
        const double sqrtw = (w[i] > 0.0 ? std::sqrt(w[i]) : 0.0);
        b[i] = sqrtw*Y[i];
        for (int j = 0; j < n_parameters; ++j)
            Acopy[i + j*n_samples] = sqrtw*A[j + i*n_parameters];
    }

    // solve

    int one(1), info;
//    const double rcond(2.2204460492503131e-16); // GSL_DBL_EPSILON
    const double rcond(1.0e-6);
    dgelsd(&n_samples, &n_parameters, &one,
           Acopy, &n_samples, b, &n_samples, S, &rcond, &rank,
           work, &lwork, iwork, &info);

    std::copy(b, b + n_parameters, x);

    chisq = 0.0;
    for (int i = 0; i < n_samples; ++i) {
        double Y_est(0);
        for (int j = 0; j < n_parameters; ++j)
            Y_est += A[j + n_parameters*i]*b[j];

        const double r = Y[i] - Y_est;
        chisq += w[i]*r*r;
    }

    // clean up

    delete[] Acopy;

    if (info != 0)
        throw std::runtime_error("MKL::dgelsd() failed");
}
#endif // USE_MKL

#ifdef USE_GSL
void wlsq::solve(int n_samples, int n_parameters,
     const double* A, // n_samples x n_parameters (row-major)
     const double* Y, // n_samples
     const double* w, // n_samples
           double* x, // n_parameters,
           double& chisq, int& rank)
{
    assert(n_samples > n_parameters);

    gsl_matrix_const_view Av =
        gsl_matrix_const_view_array(A, n_samples, n_parameters);

    gsl_vector_const_view Yv = gsl_vector_const_view_array(Y, n_samples);
    gsl_vector_const_view wv = gsl_vector_const_view_array(w, n_samples);

    gsl_vector_view xv = gsl_vector_view_array(x, n_parameters);

    // use A for the covariance matrix
    gsl_matrix* cov = gsl_matrix_alloc(n_parameters, n_parameters);

    gsl_multifit_linear_workspace* ws
        = gsl_multifit_linear_alloc(n_samples, n_parameters);

    size_t srank;
    int status = gsl_multifit_wlinear_usvd
        (&(Av.matrix), &(wv.vector), &(Yv.vector),
        GSL_DBL_EPSILON, &srank, &(xv.vector), cov, &chisq, ws);

    gsl_multifit_linear_free(ws);
    gsl_matrix_free(cov);

    rank = int(srank);
    if (status != 0)
        throw std::runtime_error("gsl_multifit_wlinear_usvd() failed");
}
#endif // USE_GSL

////////////////////////////////////////////////////////////////////////////////

} // namespace kit

////////////////////////////////////////////////////////////////////////////////
