#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[4], const double x[1],
                        double g[1])
{
    const double t1 = a[0];
    const double t2 = a[1];
    const double t6 = x[0];
    const double t4 = a[3]*t6;
    const double t5 = a[2];
    const double t7 = (t4+t5)*t6;
    const double t9 = (t2+t7)*t6;
    g[0] = ((2.0*t4+t5)*t6+t2+t7)*t6+t1+t9;
    return (t1+t9)*t6;

}

} // namespace mb_system