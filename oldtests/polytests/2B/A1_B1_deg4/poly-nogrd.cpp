#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[4], const double x[1])
{
    const double t4 = x[0];
    return((a[0]+(a[1]+(a[3]*t4+a[2])*t4)*t4)*t4);
}

} // namespace mb_system