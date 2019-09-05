#include <cmath>
#include <cassert>

#include "tang-toennies.h"

namespace x2o {

double tang_toennies(int n, const double& x)
{
    assert(n >= 0);

    int nn = n;

    double sum = 1.0 + x/nn;
    while (--nn != 0)
        sum = 1.0 + sum*x/nn;

    double tt = 1.0 - sum*std::exp(-x);

    if (std::fabs(tt) < 1.0e-8) {

        double term(1);
        for (nn = n; nn != 0; --nn)
            term *= x/nn;

        sum = 0.0;
        for (nn = n + 1; nn < 1000; ++nn) {
            term *= x/nn;
            sum += term;

            if (std::fabs(term/sum) < 1.0e-8)
                break;
        }

        tt = sum*std::exp(-x);
    }

    return tt;
}

} // namespace x2o
