#ifndef RWLSQ_H
#define RWLSQ_H

namespace kit {

//
// regularized weighted linear squares (assuming that the
// number of samples is greater than the number of unknowns)
//

struct rwlsq {
    static void solve(int n_samples, int n_parameters,
                      const double* A, // n_samples x n_parameters (row-major)
                      const double* Y,   // n_samples
                      const double* w,     // n_samples
                      const double& alpha,   // >= 0.0
                            double* x,        // n_parameters,
                            double& chisq,     // sum w_i*(Y_i - A*x)^2
                            double& penaltysq); // sum alpha^2*(x_i - <x>)^2
private:
    rwlsq();
};

} // namespace kit

#endif // RWLSQ_H
