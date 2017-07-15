#ifndef X2B_A1B2Z2_D1E2_V1_H
#define X2B_A1B2Z2_D1E2_V1_H
 
#include "poly_2b_A1B2Z2_D1E2_v1x.h" 
#include "poly_2b_A1B2Z2_D1E2.h" 
 

#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace x2b_A1B2Z2_D1E2 {

//----------------------------------------------------------------------------//

struct x2b_A1B2Z2_D1E2_v1x { 

    typedef mb_system::poly_model poly_type;
    typedef mb_system_fit::poly_model poly_type_fit;


    static std::string name();
    void load_netcdf(const char*);

    // returns 2B contribution only
    // XYZ is for the real sites
//    double operator()(const double xyz[18], double grd[18]) const; 
    double operator()(const double xyz[18]) const; 

    double eval(const double* mon1, const double* mon2) const;
    double eval(const double* mon1, const double* mon2, double* g1, double* g2) const;
    void cart_to_vars(const double xyz[18], double* vars, double& s, double& gs) const;
    
    // For fitting purposes
    size_t get_nvars() {return poly_type::n_vars;}
    size_t get_num_nonlinear_params() {return 20;}
    void set_nonlinear_parameters(const double*);
    void set(const double* xxx);
    void get_nonlinear_parameters(double*);
    bool nonlinear_parameters_out_of_range() const;
    inline double basis(const double xyz[18], double*) const; //basis polynomials
    inline size_t nparams() {
      return poly_type::size;
    }
    

    inline void as_cdl(std::ostream&) const;
    void write_cdl(std::ostream&, unsigned, unsigned, const double*) const;
    
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

// For fitting
inline void x2b_A1B2Z2_D1E2_v1x::as_cdl(std::ostream& os) const
{
    write_cdl(os, 3, poly_type::size, m_poly);
}


inline double
x2b_A1B2Z2_D1E2_v1x::basis(const double xyz[18], double* mmm) const //mmm = monomials?
{
    double v[poly_type::n_vars], s, gs;

    cart_to_vars(xyz, v, s, gs);

    poly_type_fit polyn;
    polyn.eval(v, mmm);
    for (unsigned n = 0; n < poly_type::size; ++n)
        mmm[n] *= s;

    return 0; // we will substract electrostatics and dispersion in the main
}

//----------------------------------------------------------------------------//

} // namespace x2b_A1B2Z2_D1E2

////////////////////////////////////////////////////////////////////////////////

#endif 
