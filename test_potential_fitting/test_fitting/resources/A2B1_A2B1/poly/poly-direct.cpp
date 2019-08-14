#include "poly-model.h"

namespace mb_system {

double poly_model::eval_direct(const double a[3], const double x[15])
{
    double p[3];
    p[0] = x[14];
    p[1] = x[11] + x[13] + x[10] + x[12];
    p[2] = x[6] + x[9] + x[8] + x[7];

    double energy(0);
    for(int i = 0; i < 3; ++i)
        energy += p[i]*a[i];

    return energy;

}
} // namespace mb_system