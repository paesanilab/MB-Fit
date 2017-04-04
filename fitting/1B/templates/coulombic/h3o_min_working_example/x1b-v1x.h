#ifndef X1B_V1X_H
#define X1B_V1X_H

#include "poly-1b-v1x.h"

////////////////////////////////////////////////////////////////////////////////

namespace h3o {

//----------------------------------------------------------------------------//

struct x1b_v1x {

    static const unsigned ncoeffs = poly_1b_v1x::size;

    double operator()(const double* w1) const;

    double operator()(const double* w1, double* g1) const;

    void load_poly_dat(const char*);

private:

    double m_kOH_intra;
    double m_kHH_intra;

    double m_dOH_intra;
    double m_dHH_intra;

private:
    double m_coeffs[ncoeffs];

//private:
//    double f_switch(const double& r, double& g) const;
};

//----------------------------------------------------------------------------//

} // namespace h3o

////////////////////////////////////////////////////////////////////////////////

#endif // X3B_V2X_H
