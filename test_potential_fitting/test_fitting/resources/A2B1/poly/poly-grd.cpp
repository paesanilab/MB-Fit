#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[6], const double x[3],
                        double g[3])
{
    const double t1 = a[2];
    const double t4 = x[2];
    const double t2 = t1*t4;
    const double t3 = a[0];
    const double t5 = x[1];
    const double t6 = t1*t5;
    const double t7 = a[3];
    const double t8 = t7*t4;
    const double t10 = x[0];
    const double t12 = a[4]*t10;
    const double t13 = a[5];
    const double t14 = t13*t5;
    const double t15 = t13*t4;
    const double t16 = a[1];
    const double t22 = t13*t10;
    g[0] = 2.0*t12+t14+t15+t16;
    g[1] = t22+2.0*t6+t8+t3;
    g[2] = t7*t5+2.0*t2+t22+t3;
    return (t2+t3)*t4+(t6+t8+t3)*t5+(t12+t14+t15+t16)*t10;

}

} // namespace mb_system