#include "poly-model.h"

namespace mb_system {

double poly_model::eval_direct(const double a[4], const double x[1])
{
    double p[4];
    p[0] = x[0];

    p[1] = x[0]*x[0];

    p[2] = x[0]*x[0]*x[0];

    p[3] = x[0]*x[0]*x[0]*x[0];

    double energy(0);
    for(int i = 0; i < 4; ++i)
        energy += p[i]*a[i];

    return energy;

}
} // namespace mb_system