#ifndef X2B_A1B2Z2_D1E2_V1X_H
#define X2B_A1B2Z2_D1E2_V1X_H
 
#include "poly_2b_A1B2Z2_D1E2_v1x.h" 

#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace x2b_A1B2Z2_D1E2 {

//----------------------------------------------------------------------------//

struct x2b_A1B2Z2_D1E2_v1x { 

    typedef mb_system::poly_model poly_type;


    static std::string name();
    void load_netcdf(const char*);

    // returns 2B contribution only
    // XYZ is for the real sites
    double operator()(const double xyz[18], double grd[18]) const; 
    double operator()(const double xyz[18]) const; 

    double eval(const double* mon1, const double* mon2) const;
    double eval(const double* mon1, const double* mon2, double* g1, double* g2) const;
    
private:

    double m_d_intra_AB;
    double m_k_intra_AB;
    double m_d_intra_BB;
    double m_k_intra_BB;
    double m_d_intra_DE;
    double m_k_intra_DE;
    double m_d_intra_EE;
    double m_k_intra_EE;
    double m_d_AD;
    double m_k_AD;
    double m_d_AE;
    double m_k_AE;
    double m_d_BD;
    double m_k_BD;
    double m_d_BE;
    double m_k_BE;
    double m_d_DZ;
    double m_k_DZ;
    double m_d_EZ;
    double m_k_EZ;

protected:
    double m_r2i = 6.0;
    double m_r2f = 7.0;

    double f_switch(const double&, double&) const; // At1_a -- At1_b separation

private:
    double m_poly[poly_type::size];
};

//----------------------------------------------------------------------------//

} // namespace x2b_A1B2Z2_D1E2

////////////////////////////////////////////////////////////////////////////////

#endif 
