#ifndef WLSQ_H
#define WLSQ_H

namespace kit {

//
// weighted linear squares assuming that the number
// of samples is greater than the number of unknowns
//

struct wlsq {
    static const char* implementation();

    static void solve(int n_samples, int n_parameters,
                      const double* A, // n_samples x n_parameters (row-major)
                      const double* Y, // n_samples
                      const double* w, // n_samples
                            double* x, // n_parameters,
                            double& chisq, // sum w_i*(Y_i - A*x)^2
                            int& rank);
private:
    wlsq();
};

} // namespace kit

#endif // WLSQ_H
