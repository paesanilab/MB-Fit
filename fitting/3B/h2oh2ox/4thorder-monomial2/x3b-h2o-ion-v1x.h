#ifndef X3B_H2O_ION_V1X_H
#define X3B_H2O_ION_V1X_H

#include "poly-3b-h2o-ion-v1x.h"

////////////////////////////////////////////////////////////////////////////////

namespace h2o_ion {

//----------------------------------------------------------------------------//

struct x3b_h2o_ion_v1x {

    static const unsigned ncoeffs = poly_3b_h2o_ion_v1x::size;

    double operator()(const double* w1,
                      const double* w2,
                      const double* x) const;

/*    double operator()(const double* w1,
                      const double* w2,
                      const double* x,
                      double* g1,
                      double* g2,
                      double* g3) const;
*/
    void load_netcdf(const char*);

private:
    double m_r3i;
    double m_r3f;

    double m_kOH_intra;
    double m_kHH_intra;

    double m_kOO;
    double m_kOH;
    double m_kHH;

    double m_kXO;
    double m_kXH;

    double m_dOH_intra;
    double m_dHH_intra;

    double m_dOO;
    double m_dOH;
    double m_dHH;

    double m_dXO;
    double m_dXH;

private:
    double m_coeffs[ncoeffs];

private:
    double f_switch(const double& r, double& g) const;
};

//----------------------------------------------------------------------------//

} // namespace h2o_ion

////////////////////////////////////////////////////////////////////////////////

#endif // X3B_H2O_ION_V1X_H
