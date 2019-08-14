#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[6], const double x[3])
{
    const double t1 = a[2];
    const double t3 = a[0];
    const double t13 = a[5];
    const double t2 = x[2];
    const double t7 = x[1];
    const double t17 = x[0];
    return((t1*t2+t3)*t2+(t1*t7+a[3]*t2+t3)*t7+(t13*t2+t13*t7+a[4]*t17+a[1])*t17);

}

} // namespace mb_system