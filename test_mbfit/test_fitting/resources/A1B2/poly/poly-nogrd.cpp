#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[6], const double x[3])
{
    const double t6 = a[3];
    const double t2 = x[2];
    const double t9 = a[2]*t2;
    const double t10 = a[0];
    const double t11 = x[1];
    const double t15 = x[0];
    return((a[5]*t2+a[1])*t2+(t6*t11+t10+t9)*t11+(a[4]*t11+t6*t15+t10+t9)*t15);
}

} // namespace mb_system