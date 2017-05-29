
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <netcdf.h>

#include "stuff.h"
#include "constants.h"
#include "x2b_A1B2_A1B2.h"

#ifdef DEBUG
#define PR(x) std::cout << #x << ": " << (x) << std::endl;
#else
#define PR(x)
#endif /* DEBUG */

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int kode)
{
    std::cerr << " ** Fatal Error in x2b_v9<D>::load_netcdf() ** : "
              << nc_strerror(kode) << std::endl;
    std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------//

void xpoints(const double& in_plane_g,
             const double& out_of_plane_g,
             const double xyz[9],
             double xa[3], double xb[3])   //calculates lone pair sites coordinates
{
    double e1[3], e2[3];

    for (int i = 0; i < 3; ++i) {
        e1[i] = xyz[i + 3] - xyz[i];
        e2[i] = xyz[i + 6] - xyz[i];
    }

    const double v1[3] = {
        (e1[0] + e2[0])/2,
        (e1[1] + e2[1])/2,
        (e1[2] + e2[2])/2
    };

    const double v2[3] = {
        e1[1]*e2[2] - e1[2]*e2[1],
        e1[2]*e2[0] - e1[0]*e2[2],
        e1[0]*e2[1] - e1[1]*e2[0]
    };

    for (int i = 0; i < 3; ++i) {
        const double in_plane = xyz[i] + in_plane_g*v1[i];
        const double out_of_plane = out_of_plane_g*v2[i];

        xa[i] = in_plane + out_of_plane;
        xb[i] = in_plane - out_of_plane;
    }
}

//----------------------------------------------------------------------------//

inline double var_intra(const double& d0, const double& k, const double& d)
{
    return std::exp(k*(d0 - d));
}

inline double var_coul(const double& d0, const double& k, const double& d)
{
    return std::exp(k*(d0 - d))/d;
}

// var_lp?
inline double var_main(const double& d0, const double& k, const double& d)
{
    return std::exp(k*(d0 - d));
}

//----------------------------------------------------------------------------//

const double delta_scale =  1.0;

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace tset {

//----------------------------------------------------------------------------//

double x2b_A1B2_A1B2::f_switch(const double& r) const 

{
    if (r > m_r2f) {
        return 0.0;
    } else if (r > m_r2i) {
        const double x = (r - m_r2i)/(m_r2f - m_r2i);
        return (1.0 + std::cos(M_PI*x))/2;
    } else {
        return 1.0;
    }
}

//----------------------------------------------------------------------------//
void x2b_A1B2_A1B2::set_nonlinear_parameters(const double* xxx) 
 { 
     m_k_intra_AB = *xxx++; 
    m_d_intra_AB = *xxx++; 
    m_k_intra_BB = *xxx++; 
    m_d_intra_BB = *xxx++; 
    m_k_AA = *xxx++; 
    m_d_AA = *xxx++; 
    m_k_AB = *xxx++; 
    m_d_AB = *xxx++; 
    m_k_BB = *xxx++; 
    m_d_BB = *xxx++; 
} 
void x2b_A1B2_A1B2::get_nonlinear_parameters(double* xxx) const { 
    PR(m_k_intra_AB) 
 
    PR(m_d_intra_AB) 
 
    PR(m_k_intra_BB) 
 
    PR(m_d_intra_BB) 
 
    PR(m_k_AA) 
 
    PR(m_d_AA) 
 
    PR(m_k_AB) 
 
    PR(m_d_AB) 
 
    PR(m_k_BB) 
 
    PR(m_d_BB) 
 
    *xxx++ = m_k_intra_AB; 
    *xxx++ = m_d_intra_AB; 
    *xxx++ = m_k_intra_BB; 
    *xxx++ = m_d_intra_BB; 
    *xxx++ = m_k_AA; 
    *xxx++ = m_d_AA; 
    *xxx++ = m_k_AB; 
    *xxx++ = m_d_AB; 
    *xxx++ = m_k_BB; 
    *xxx++ = m_d_BB; 

}

//----------------------------------------------------------------------------//

bool x2b_A1B2_A1B2::nonlinear_parameters_out_of_range() const { 

    // ##DEFINE HERE## the maximum and minimum for the ks and ds
    const double k_min =  0.0;
    const double k_max = 3.0;
    const double k_min_intra =  0.0;
    const double k_max_intra= 2.0;

    const double d_min = 0.0;
    const double d_max =  7.0;
    const double d_min_intra =  0.0;
    const double d_max_intra= 3.0;
    
return false
       || m_k_intra_AB < k_min_intra 
       || m_k_intra_AB > k_max_intra 
       || m_d_intra_AB < d_min_intra 
       || m_d_intra_AB > d_max_intra 
       || m_k_intra_BB < k_min_intra 
       || m_k_intra_BB > k_max_intra 
       || m_d_intra_BB < d_min_intra 
       || m_d_intra_BB > d_max_intra 
       || m_k_AA < k_min 
       || m_k_AA > k_max 
       || m_d_AA < d_min 
       || m_d_AA > d_max 
       || m_k_AB < k_min 
       || m_k_AB > k_max 
       || m_d_AB < d_min 
       || m_d_AB > d_max 
       || m_k_BB < k_min 
       || m_k_BB > k_max 
       || m_d_BB < d_min 
       || m_d_BB > d_max ; 


}

//----------------------------------------------------------------------------//


void x2b_A1B2_A1B2::do_setup(const std::map<std::string, double>& kv) { 
    m_k_intra_AB = 0.5;
    m_d_intra_AB = 0.5;
    m_k_intra_BB = 0.5;
    m_d_intra_BB = 0.5;
    m_k_AA = 0.5;
    m_d_AA = 0.5;
    m_k_AB = 0.5;
    m_d_AB = 0.5;
    m_k_BB = 0.5;
    m_d_BB = 0.5;
    m_r2i = 6.0;
    m_r2f = 7.0;


    // override from the KV

    std::map<std::string, double>::const_iterator iter;

#   define OVERRIDE(name)       \
    iter = kv.find(#name);      \
    if (iter != kv.end())       \
        m_##name = iter->second;

    OVERRIDE(k_intra_AB);
    OVERRIDE(d_intra_AB);
    OVERRIDE(k_intra_BB);
    OVERRIDE(d_intra_BB);
    OVERRIDE(k_AA);
    OVERRIDE(d_AA);
    OVERRIDE(k_AB);
    OVERRIDE(d_AB);
    OVERRIDE(k_BB);
    OVERRIDE(d_BB);

#   undef OVERRIDE
}

//----------------------------------------------------------------------------//

// ##DEFINE HERE## need to think how to generalize this function...

void  x2b_A1B2_A1B2::cart_to_vars(const double* xyz, double* v, double& s) const { 

    const double* A = xyz;
    const double* B1 = xyz + 3;
    const double* B2 = xyz + 6;
    
    const double* A_2 = xyz + 9;
    const double* B1_2 = xyz + 12;
    const double* B2_2 = xyz + 15;
    
    
// ##DEFINE HERE## the lone pairs if any
    // double Xa1[3], Xa2[3];

    // xpoints(m_in_plane_gamma, m_out_of_plane_gamma, O, Xa1, Xa2);


    using x2o::distance;

    v[0]  = var_intra(m_d_intra_BB, m_k_intra_BB, distance(B1, B2));
    v[1]  = var_intra(m_d_intra_BB, m_k_intra_BB, distance(B1_2, B2_2));
    v[2]  = var_intra(m_d_intra_AB, m_k_intra_AB, distance(A, B1));
    v[3]  = var_intra(m_d_intra_AB, m_k_intra_AB, distance(A, B2));
    v[4]  = var_intra(m_d_intra_AB, m_k_intra_AB, distance(A_2, B1_2));
    v[5]  = var_intra(m_d_intra_AB, m_k_intra_AB, distance(A_2, B2_2));
    
    v[6]  = var_intra(m_d_BB, m_k_BB, distance(B1, B1_2));
    v[7]  = var_intra(m_d_BB, m_k_BB, distance(B1, B2_2));
    v[8]  = var_intra(m_d_BB, m_k_BB, distance(B2, B1_2));
    v[9]  = var_intra(m_d_BB, m_k_BB, distance(B2, B2_2));
    
    v[10]  = var_intra(m_d_AB, m_k_AB, distance(A, B1_2));
    v[11]  = var_intra(m_d_AB, m_k_AB, distance(A, B2_2));
    v[12]  = var_intra(m_d_AB, m_k_AB, distance(A_2, B1));
    v[13]  = var_intra(m_d_AB, m_k_AB, distance(A_2, B2));
    
    v[14]  = var_intra(m_d_AA, m_k_AA, distance(A, A_2));

    s = f_switch(distance(A, A_2));

  PR(s)
  PR(v[0]);
  PR(v[1]);
  PR(v[2]);
  PR(v[3]);
  PR(v[4]);
  PR(v[5]);
  PR(v[6]);
  PR(v[7]);
  PR(v[8]);
  PR(v[9]);
  PR(v[10]);
  PR(v[11]);
  PR(v[12]);
  PR(v[13]);
  PR(v[14]);
} 

//----------------------------------------------------------------------------//

void x2b_A1B2_A1B2::write_cdl(std::ostream& os, unsigned DEG,
                                   unsigned npoly, const double* poly) const {

    ios::fmtflags saved_flags = os.flags();

    os << setprecision(15) << scientific
       << "netcdf x2b_A1B2_A1B2{ " << endl 
       << "  // global attributes" << endl 
       << "  :name = x2b_A1B2_A1B2<" << DEG << ">;" << endl; 

     x2b_A1B2_A1B2::as_cdl(os); 
 
    os        << "  :k_intra_AB = " << setw(22) << m_k_intra_AB << "; // A^(-1))" << endl 
       << "  :d_intra_AB = " << setw(22) << m_d_intra_AB << "; // A^(-1))" << endl 
       << "  :k_intra_BB = " << setw(22) << m_k_intra_BB << "; // A^(-1))" << endl 
       << "  :d_intra_BB = " << setw(22) << m_d_intra_BB << "; // A^(-1))" << endl 
       << "  :k_AA = " << setw(22) << m_k_AA << "; // A^(-1))" << endl 
       << "  :d_AA = " << setw(22) << m_d_AA << "; // A^(-1))" << endl 
       << "  :k_AB = " << setw(22) << m_k_AB << "; // A^(-1))" << endl 
       << "  :d_AB = " << setw(22) << m_d_AB << "; // A^(-1))" << endl 
       << "  :k_BB = " << setw(22) << m_k_BB << "; // A^(-1))" << endl 
       << "  :d_BB = " << setw(22) << m_d_BB << "; // A^(-1))" << endl 

          "  :r2i = " << setw(22) << m_r2i << "; // A" << endl
          "  :r2f = " << setw(22) << m_r2f << "; // A" << endl
          "  dimensions:" << endl
          "  poly = " << npoly << ";" << endl
          "  variables:" << endl
          "    double poly(poly);" << endl
          "data:" << endl;

    os << "poly =" << endl;

    {
        unsigned n(0);
        for (n = 0; n + 1 < npoly; ++n)
            os << std::setw(22) << poly[n] << ", // " << n << endl;
        os <<  std::setw(22) << poly[n] << "; // " << n << endl << endl;
    }

    os << "}" << end;
    os.flags(saved_flags);
} 

//----------------------------------------------------------------------------//

void x2b_A1B2_A1B2::do_load_netcdf 
    (const char* fn, unsigned npoly, double* poly) { 
}
    assert(fn);

    int rc, ncid;
    if ((rc = nc_open(fn, NC_NOWRITE, &ncid)))
        error(rc);

    from_cdf(ncid);

#   define RETRIEVE(name) \
    if ((rc = nc_get_att_double(ncid, NC_GLOBAL, #name , &m_##name))) \
        error(rc);
        
RETRIEVE(k_intra_AB)
RETRIEVE(d_intra_AB)
RETRIEVE(k_intra_BB)
RETRIEVE(d_intra_BB)
RETRIEVE(k_AA)
RETRIEVE(d_AA)
RETRIEVE(k_AB)
RETRIEVE(d_AB)
RETRIEVE(k_BB)
RETRIEVE(d_BB)


#   undef RETRIEVE

    int varid;

    if ((rc = nc_inq_varid(ncid, "poly", &varid)))
        error(rc);

    for (size_t n = 0; n < npoly; ++n) {
        if ((rc = nc_get_var1_double(ncid, varid, &n, poly + n)))
            error(rc);
    }

    if ((rc = nc_close(ncid)))
        error(rc);
}

//----------------------------------------------------------------------------//

} // namespace tset

////////////////////////////////////////////////////////////////////////////////
