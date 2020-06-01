#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[6], const double x[3],
                        double g[3])
{
    const double t4 = x[2];
    const double t2 = a[5]*t4;
    const double t3 = a[1];
    const double t6 = a[3];
    const double t5 = x[1];
    const double t7 = t6*t5;
    const double t8 = a[2];
    const double t9 = t8*t4;
    const double t10 = a[0];
    const double t11 = x[0];
    const double t13 = t6*t11;
    const double t14 = a[4];
    const double t15 = t14*t5;
    g[0] = 2.0*t13+t15+t9+t10;
    g[1] = t11*t14+t10+2.0*t7+t9;
    g[2] = t11*t8+t5*t8+2.0*t2+t3;
    return (t2+t3)*t4+(t7+t9+t10)*t5+(t13+t15+t9+t10)*t11;

}

} // namespace mb_system