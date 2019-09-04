#ifdef HAVE_CONFIG_H
#   include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iostream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

//#undef USE_MKL
#define USE_MKL 1 // FIXME remove this 
#ifdef USE_MKL
#   include <mkl_lapack.h>
#endif // USE_MKL

#include "rwlsq.h"

#ifdef USE_MKL

namespace {

// dgesdd(const char* jobz, const MKL_INT* m, const MKL_INT* n, double* a,
//        const MKL_INT* lda, double* s, double* u, const MKL_INT* ldu,
//        double* vt, const MKL_INT* ldvt, double* work,
//        const MKL_INT* lwork, MKL_INT* iwork, MKL_INT* info);

void query_workspace_sizes(int M, int N, int& lwork, int& liwork)
{
    const char jobz = 'O';
    int minus_one(-1), info;
    double work;

    lwork = -1;
    dgesdd(&jobz, &M, &N,
           0, &M, // A
           0,     // S
           0, &M, // U
           0, &N, // Vt
           &work, &minus_one,
           0, &info);

    assert(info == 0);
    lwork = int(work);
    liwork = 8*std::min(M, N);
}

void invoke_dgesdd(gsl_matrix* A, gsl_matrix* V, gsl_vector* S)
{
    const int M = A->size1;
    const int N = A->size2;

    int lwork, liwork;
    query_workspace_sizes(M, N, lwork, liwork);

    const size_t nmem = M*N
                      + lwork
                      + liwork*sizeof(double)/sizeof(int);

    double* mem = new double[nmem];

    double* A_col_major = mem;
    double* work = A_col_major + M*N;
    int* iwork = reinterpret_cast<int*>(work + lwork);

    // populate A_col_major

    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            A_col_major[j*M + i] = gsl_matrix_get(A, i, j);

    const char jobz = 'O';
    int info;

    dgesdd(&jobz, &M, &N,
           A_col_major, &M, // A
           S->data,     // S
           0, &M,       // U
           V->data, &N, // Vt
           work, &lwork,
           iwork, &info);

    // transpose A_col_major into A

    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            gsl_matrix_set(A, i, j, A_col_major[j*M + i]);

    delete[] mem;

    if (info != 0)
        throw std::runtime_error("MKL::dgesdd() failed");
}

} // namespace

#endif // USE_MKL

namespace kit {

void rwlsq::solve(int n_samples, int n_parameters,
                  const double* A,   // n_samples x n_parameters (row-major)
                  const double* Y,    // n_samples
                  const double* w,     // n_samples
                  const double& alpha,  // > 0.0
                        double* x,       // n_parameters,
                        double& chisq,    // sum w_i*(Y_i - A*x)^2
                        double& penaltysq) // sum alpha^2*(x_i - <x>)^2 ??
{
//    std::cout << "Entered rwlsq::solve() \n" ;

    assert(n_samples > n_parameters);
    assert(alpha > 0.0);

    gsl_matrix_const_view A_view =
        gsl_matrix_const_view_array(A, n_samples, n_parameters);

    // scale A: As = sqrt(w) A

    gsl_matrix* As = gsl_matrix_alloc(n_samples, n_parameters);

    gsl_matrix_memcpy(As, &A_view.matrix);

    for (int i = 0; i < n_samples; ++i) {
        double wi = w[i];

        if (wi < 0.0)
            wi = 0.0;

        {
            gsl_vector_view row = gsl_matrix_row(As, i);
            gsl_vector_scale(&row.vector, std::sqrt(wi));
        }
    }

    // decompose As into U S Q^T

    gsl_matrix* Q = gsl_matrix_alloc(n_parameters, n_parameters);
    gsl_matrix* QSI = gsl_matrix_alloc(n_parameters, n_parameters);

    gsl_vector* S = gsl_vector_alloc(n_parameters);
    gsl_vector* xt = gsl_vector_alloc(n_parameters); // FIXME why was it _calloc?

#   ifdef USE_MKL
    invoke_dgesdd(As, Q, S);
#   else
//    std::cout << "Going to enter gsl_linalg_SV_decomp_mod()" << std::endl;
    gsl_linalg_SV_decomp_mod(As, QSI, Q, S, xt);
#   endif // USE_MKL

    // solve sqrt(w) Y = As x for x, by first computing t = sqrt(w) Y */

    gsl_vector* t = gsl_vector_alloc(n_samples);

    for (int i = 0; i < n_samples; ++i) {
        double wi = w[i];
        double yi = Y[i];

        if (wi < 0.0)
            wi = 0.0;

        gsl_vector_set(t, i, std::sqrt(wi)*yi);
    }

    gsl_blas_dgemv(CblasTrans, 1.0, As, t, 0.0, xt);
    gsl_vector_free(t);
    gsl_matrix_free(As);

    // scale the matrix Q,  Q' = Q S^-1

    gsl_matrix_memcpy(QSI, Q);
    gsl_matrix_free(Q);

    const double alpha2 = alpha*alpha;
    for (int j = 0; j < n_parameters; ++j) {
        gsl_vector_view column = gsl_matrix_column(QSI, j);
        const double sigma = gsl_vector_get(S, j);

        const double fj = sigma/(sigma*sigma + alpha2);

        gsl_vector_scale(&column.vector, fj);
    }

    std::fill(x, x + n_parameters, 0.0);
    gsl_vector_view x_view = gsl_vector_view_array(x, n_parameters);

    // solution

    gsl_blas_dgemv(CblasNoTrans, 1.0, QSI, xt, 0.0, &x_view.vector);
    gsl_vector_free(xt);
    gsl_vector_free(S);
    gsl_matrix_free(QSI);

    // compute chisq from residual r = Y - A x

    chisq = 0.0;

    for (int i = 0; i < n_samples; ++i) {
        const double yi = Y[i];
        double wi = w[i];

        if (wi < 0.0)
            wi = 0.0;

        gsl_vector_const_view row =
            gsl_matrix_const_row(&A_view.matrix, i);

        double y_est;
        gsl_blas_ddot(&row.vector, &x_view.vector, &y_est);

        const double ri = yi - y_est;
        chisq += wi*ri*ri;
    }

    penaltysq = gsl_blas_dnrm2(&x_view.vector); //\sqrt {\sum x_i^2}
    penaltysq = alpha2*penaltysq*penaltysq;
}

} // namespace kit
