#ifndef X2B_H2O_ION_V1X_H
#define X2B_H2O_ION_V1X_H

#include "x2b-h2o-ion-base.h"
#include "poly-2b-h2o-ion-v1x.h"

////////////////////////////////////////////////////////////////////////////////

namespace h2o_ion {

//----------------------------------------------------------------------------//

//
//  x2b_h2o_ion_v1<4> with gradients and no fitting interface
//

//----------------------------------------------------------------------------//

struct x2b_h2o_ion_v1x : public x2b_h2o_ion_base {
    typedef x2o::poly_2b_h2o_ion_v1x poly_type;

    static std::string name();
    void load_netcdf(const char*);

    // returns 2B contribution only
    double operator()(const double xyz[12], double grd[12]) const;

    double polynomial(const double* w1, const double* x) const;
//    double polynomial(const double* w1, const double* x, double* g1, double* g2) const;

private:
    double m_k_HH_intra;
    double m_k_OH_intra;

    double m_k_XH_coul;
    double m_k_XO_coul;

    double m_k_xlp_main;

    double m_in_plane_gamma;
    double m_out_of_plane_gamma;

protected:
    double m_r2i;
    double m_r2f;

    double f_switch(const double&, double&) const; // O-X separation

private:
    double m_poly[poly_type::size];
};

//----------------------------------------------------------------------------//

} // namespace h2o_ion

////////////////////////////////////////////////////////////////////////////////

#endif // X2B_H2O_ION_V1X_H
