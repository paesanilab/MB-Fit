#ifndef X2B_A1B2_A1B2_H
#define X2B_A1B2_A1B2_H

#include "poly_2b_A1B2_A1B2.h"

#include <map>
#include <string>
#include <iostream>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////

namespace tset {

//----------------------------------------------------------------------------//
template <unsigned DEG> struct x2b { //template class based on degree of the polynomial
  public:

typedef tset::poly_2b_A1B2_A1B2<DEG> poly_type;

    static inline std::string name();

    // fitting interface
    inline static size_t nparams(); // number of fitting parameters
    inline double basis(const double xyz[18], double*) const; //basis polynomials


    // value
    inline void set(const double*);
    inline double value(const double xyz[18]) const;


    // IO
    inline void as_cdl(std::ostream&) const;
    inline void load_netcdf(const char*);
    
    x2b()
    {
        poly = new double[poly_type::size];
    }

    ~x2b()
    {
        delete[] poly;
    }
    
    static const unsigned num_nonlinear_params = 10;

    void set_nonlinear_parameters(const double*); //modify these functions in cpp
    void get_nonlinear_parameters(double*) const;

    bool nonlinear_parameters_out_of_range() const;
    
  protected:
    double m_k_intra_AB;
    double m_d_intra_AB;
    double m_k_intra_BB;
    double m_d_intra_BB;
    double m_k_AA;
    double m_d_AA;
    double m_k_AB;
    double m_d_AB;
    double m_k_BB;
    double m_d_BB;
    void do_setup(const std::map<std::string, double>&); //gives default values to the non linear params and the switching function params 
    void cart_to_vars(const double xyz[18], double* vars, double& s) const; //cartesian to the required distance variables 


  protected:
    double m_r2i;  //switching function parameters, inner and outer radii 
    double m_r2f;

    double f_switch(const double&) const; // 1st atom - 1st atom separation

  protected:
    void write_cdl(std::ostream&, unsigned, unsigned, const double*) const;
    void do_load_netcdf(const char* fn, unsigned, double*);

  private:
    double* poly;
};

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline size_t x2b<DEG>::nparams()
{
    return poly_type::size;
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline double
x2b<DEG>::basis(const double xyz[18], double* mmm) const 

{
    double v[poly_type::n_vars], s;

    cart_to_vars(xyz, v, s);

    poly_type::eval(v, mmm);
    for (unsigned n = 0; n < poly_type::size; ++n)
        mmm[n] *= s;

    return x2b::value(xyz);
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline void x2b<DEG>::set(const double* p)
{
    std::copy(p, p + poly_type::size, poly);
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline double x2b<DEG>::value(const double xyz[18]) const 

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

    return s*E_poly;

}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline void x2b<DEG>::as_cdl(std::ostream& os) const
{
    write_cdl(os, DEG, poly_type::size, poly);
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline void x2b<DEG>::load_netcdf(const char* fn)
{
    do_load_netcdf(fn, poly_type::size, poly);
}

//----------------------------------------------------------------------------//

template <unsigned DEG>
inline std::string x2b<DEG>::name()
{
    std::ostringstream oss;
    oss << "x2b<" << DEG << ">";
    return oss.str();
}

//----------------------------------------------------------------------------//

} // namespace tset

////////////////////////////////////////////////////////////////////////////////

#endif 
