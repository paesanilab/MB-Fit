#include "poly-model.h"

namespace mb_system {

double poly_model::eval_direct(const double a[6], const double x[3])
{
    double p[6];
    p[0] = x[1] + x[2];
    p[1] = x[0];

    p[2] = x[1]*x[1] + x[2]*x[2];
    p[3] = x[1]*x[2];
    p[4] = x[0]*x[0];
    p[5] = x[0]*x[1] + x[0]*x[2];

    double energy(0);
    for(int i = 0; i < 6; ++i)
        energy += p[i]*a[i];

    return energy;

}
} // namespace mb_system