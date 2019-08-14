#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[3], const double x[15],
                        double g[15])
{
    const double t1 = a[2];
    const double t6 = a[1];
    const double t11 = a[0];
    g[0] = 0.0;
    g[1] = 0.0;
    g[2] = 0.0;
    g[3] = 0.0;
    g[4] = 0.0;
    g[5] = 0.0;
    g[6] = t1;
    g[7] = t1;
    g[8] = t1;
    g[9] = t1;
    g[10] = t6;
    g[11] = t6;
    g[12] = t6;
    g[13] = t6;
    g[14] = t11;
    return t1*x[6]+t1*x[7]+t1*x[8]+t1*x[9]+t11*x[14]+t6*x[10]+t6*x[11]+
t6*x[12]+t6*x[13];

}

} // namespace mb_system