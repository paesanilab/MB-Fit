#include "poly-model.h"

namespace mb_system {

double poly_model::eval(const double a[3], const double x[15])
{
    const double t1 = a[2];
    const double t6 = a[1];
    return(t1*x[6]+t1*x[7]+t1*x[8]+t1*x[9]+t6*x[10]+t6*x[11]+t6*x[12]+t6*x[13]+a[0]*x[14]);

}

} // namespace mb_system