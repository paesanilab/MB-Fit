#ifndef X2B_H2O_X_V1_H
#define X2B_H2O_X_V1_H

#include <sstream>

#include "x2b-h2o-ion-base.h"   //base header

#include "poly-2b-h2o-ion-v1.h"    //polynomial

//#define DEBUG_ENERGY
//#define NO_POLY

#ifdef DEBUG_ENERGY 
#include <iostream>
#endif 
////////////////////////////////////////////////////////////////////////////////

namespace h2o_ion {

//----------------------------------------------------------------------------//

//
//  something on top of h2o_ion::x2b_h2o_ion_base
//

//----------------------------------------------------------------------------//

struct x2b_h2o_ion_v1_common : public x2b_h2o_ion_base {  //inherits stuff from base 
public:
//    static const unsigned num_nonlinear_params = 2 // gamma
//                                               + 2 //kcoul
//                                               + 1 //kmain
//					       + 2 //kintra
//					       + 2 //dcoul
//					       + 1 //dmain
//					       + 2; //dintra
    static const unsigned num_nonlinear_params = 2 //kcoul
                                               + 1 //kmain
					       + 2 //kintra
					       + 2 //dcoul
					       + 1 //dmain
					       + 2; //dintra
    void set_nonlinear_parameters(const double*); //modify these functions in cpp
    void get_nonlinear_parameters(double*) const;

    bool nonlinear_parameters_out_of_range() const;

protected:
    x2b_h2o_ion_v1_common();

protected:
    double m_k_HH_intra;
    double m_k_OH_intra;

    double m_k_XH_coul;
    double m_k_XO_coul;

    double m_k_xlp_main;

    double m_d_HH_intra;
    double m_d_OH_intra;

    double m_d_XH_coul;
    double m_d_XO_coul;

    double m_d_xlp_main;

protected:
    double m_in_plane_gamma;
    double m_out_of_plane_gamma;

    void do_setup(const std::map<std::string, double>&); //gives default values to the non linear params and the switching function params 
    void cart_to_vars(const double xyz[12], double* vars, double& s) const; //cartesian to the required distance variables

protected:
    double m_r2i;  //switching function parameters, inner and outer radii 
    double m_r2f;

    double f_switch(const double&) const; // O-O separation

protected:
    void write_cdl(std::ostream&, unsigned, unsigned, const double*) const;
    void do_load_netcdf(const char* fn, unsigned, double*);
};

//----------------------------------------------------------------------------//

template <unsigned DEG> //template class based on degree of the polynomial
struct x2b_h2o_ion_v1 : public x2b_h2o_ion_v1_common {

    typedef h2o_ion::poly_2b_h2o_ion<DEG> poly_type;

    static inline std::string name();

    // fitting interface
    inline static size_t nparams(); // number of fitting parameters
    inline double basis(const double xyz[12], double*) const; //basis polynomials

    // value
    inline void set(const double*);
    inline double value(const double xyz[12]) const;

    // IO
    inline void as_cdl(std::ostream&) const;
    inline void load_netcdf(const char*);

    x2b_h2o_ion_v1()
    {
        poly = new double[poly_type::size];
    }

    ~x2b_h2o_ion_v1()
    {
        delete[] poly;
    }

private:
    double* poly;
};

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline size_t x2b_h2o_ion_v1<DEG>::nparams()
{
    return poly_type::size;
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline double
x2b_h2o_ion_v1<DEG>::basis(const double xyz[12], double* mmm) const //mmm = monomials?
{
    double v[poly_type::n_vars], s;

    cart_to_vars(xyz, v, s);

    poly_type::eval(v, mmm);
    for (unsigned n = 0; n < poly_type::size; ++n)
        mmm[n] *= s;

    return x2b_h2o_ion_base::value(xyz);
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline void x2b_h2o_ion_v1<DEG>::set(const double* p)
{
    std::copy(p, p + poly_type::size, poly);
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline double x2b_h2o_ion_v1<DEG>::value(const double xyz[12]) const
{
    double v[poly_type::n_vars], s;
    cart_to_vars(xyz, v, s);

    double E_poly(0);

    {
        double mono[poly_type::size];
        poly_type::eval(v, mono);
		
        for (size_t n = 0; n < poly_type::size; ++n)
            E_poly += poly[n]*mono[n];
    }

#ifdef DEBUG_ENERGY 
    std::cout << "s*E_poly = " << s*E_poly << std::endl;
#endif 

#ifdef NO_POLY
    return x2b_h2o_ion_base::value(xyz);
#else
    return s*E_poly + x2b_h2o_ion_base::value(xyz);
#endif 
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline void x2b_h2o_ion_v1<DEG>::as_cdl(std::ostream& os) const
{
    write_cdl(os, DEG, poly_type::size, poly);
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline void x2b_h2o_ion_v1<DEG>::load_netcdf(const char* fn)
{
    do_load_netcdf(fn, poly_type::size, poly);
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline std::string x2b_h2o_ion_v1<DEG>::name()
{
    std::ostringstream oss;
    oss << "x2b_h2o_ion_v1<" << DEG << ">";
    return oss.str();
}

//----------------------------------------------------------------------------//

} // namespace h2o_ion

////////////////////////////////////////////////////////////////////////////////

#endif // X2B_H2O_X_V1_H
