#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <netcdf.h>

#include "stuff.h" 
#include "x2b-h2o-ion-v1.h"
#include "constants.h"

//#define DEBUG 

//#define CHLORIDE 
//#define FLUORIDE 
//#define BROMIDE  
//#define IODIDE 
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

namespace h2o_ion {

//----------------------------------------------------------------------------//

x2b_h2o_ion_v1_common::x2b_h2o_ion_v1_common()
{
    setup("");
}

//----------------------------------------------------------------------------//

double x2b_h2o_ion_v1_common::f_switch(const double& r) const //if at all, how should this be modified?
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

void x2b_h2o_ion_v1_common::set_nonlinear_parameters(const double* xxx)
{
    m_k_HH_intra = *xxx++;
    m_k_OH_intra = *xxx++;

    m_k_XH_coul = *xxx++;
    m_k_XO_coul = *xxx++;

    m_k_xlp_main = *xxx++;

    m_d_HH_intra = *xxx++;
    m_d_OH_intra = *xxx++;

    m_d_XH_coul = *xxx++;
    m_d_XO_coul = *xxx++;

    m_d_xlp_main = *xxx++;

//    m_in_plane_gamma = *xxx++;
//    m_out_of_plane_gamma = *xxx++;
}

//----------------------------------------------------------------------------//

void x2b_h2o_ion_v1_common::get_nonlinear_parameters(double* xxx) const
{
    PR(m_k_HH_intra) 
    PR(m_k_OH_intra) 
    PR(m_k_XH_coul) 
    PR(m_k_XO_coul) 
    PR(m_k_xlp_main)
    PR(m_d_HH_intra) 
    PR(m_d_OH_intra) 
    PR(m_d_XH_coul) 
    PR(m_d_XO_coul) 
    PR(m_d_xlp_main)
//    PR(m_in_plane_gamma)
//    PR(m_out_of_plane_gamma)

    *xxx++ = m_k_HH_intra;
    *xxx++ = m_k_OH_intra;

    *xxx++ = m_k_XH_coul;
    *xxx++ = m_k_XO_coul;

    *xxx++ = m_k_xlp_main;

    *xxx++ = m_d_HH_intra;
    *xxx++ = m_d_OH_intra;

    *xxx++ = m_d_XH_coul;
    *xxx++ = m_d_XO_coul;

    *xxx++ = m_d_xlp_main;

//    *xxx++ = m_in_plane_gamma;
//    *xxx++ = m_out_of_plane_gamma;
}

//----------------------------------------------------------------------------//

bool x2b_h2o_ion_v1_common::nonlinear_parameters_out_of_range() const
{
    const double k_min =  0.0;
    const double k_max = 3.0;
    const double k_min_intra =  0.0;
    const double k_max_intra= 2.0;

    const double d_min = 0.0;
    const double d_max =  7.0;
    const double d_min_intra =  0.0;
    const double d_max_intra= 3.0;

//    const double oop_gamma_max = 2.0;
//    const double ip_gamma_max = 2.0;
//
//    const double oop_gamma_min = 0.0;
//    const double ip_gamma_min = 0.0;

    // TODO: Put in reasonable bounds for nonlinear parameters
        // mins
    return m_k_HH_intra < k_min_intra || m_k_OH_intra < k_min_intra
        || m_k_XH_coul < k_min || m_k_XO_coul < k_min 
        || m_k_xlp_main < k_min || m_d_HH_intra < d_min_intra
        || m_d_OH_intra < d_min_intra || m_d_XH_coul < d_min
        || m_d_XO_coul < d_min || m_d_xlp_main < d_min
	// maxes
	|| m_k_HH_intra > k_max_intra || m_k_OH_intra > k_max_intra
        || m_k_XH_coul > k_max || m_k_XO_coul > k_max 
        || m_k_xlp_main > k_max || m_d_HH_intra > d_max_intra
        || m_d_OH_intra > d_max_intra || m_d_XH_coul > d_max
        || m_d_XO_coul > d_max || m_d_xlp_main > d_max;
	// lone-pair sites
//        || std::fabs(m_out_of_plane_gamma) > oop_gamma_max
//        || std::fabs(m_in_plane_gamma) > ip_gamma_max
//        || std::fabs(m_out_of_plane_gamma) < oop_gamma_min
//        || std::fabs(m_in_plane_gamma) < ip_gamma_min; //TODO check!
}

//----------------------------------------------------------------------------//

void x2b_h2o_ion_v1_common::do_setup(const std::map<std::string, double>& kv)
{
    x2b_h2o_ion_base::do_setup(kv);

    // set defaults

    m_k_HH_intra = 0.4;
    m_k_OH_intra = 0.6;

    m_k_XH_coul = 0.8; //still using the old values for water dimer
    m_k_XO_coul = 0.8;

    m_k_xlp_main = 0.6; //old values for water dimer 

    m_d_HH_intra = 0.6;
    m_d_OH_intra = 0.4;

    m_d_XH_coul = 2.0; 
    m_d_XO_coul = 2.0;

    m_d_xlp_main = 2.0;

    m_in_plane_gamma =     -9.721486914088159e-02;  //from MBpol
    m_out_of_plane_gamma = 9.859272078406150e-02;

#ifdef FLUORIDE
    m_r2i = 5.0; 
    m_r2f = 6.0;
#endif 
#ifdef CHLORIDE
    m_r2i = 5.5; 
    m_r2f = 6.5;
#endif 
#ifdef BROMIDE 
    m_r2i = 5.5; 
    m_r2f = 6.5;
#endif 
#ifdef IODIDE  
    m_r2i = 6.0; 
    m_r2f = 7.0;
#endif 
    // override from the KV

    std::map<std::string, double>::const_iterator iter;

#   define OVERRIDE(name)       \
    iter = kv.find(#name);      \
    if (iter != kv.end())       \
        m_##name = iter->second;

    OVERRIDE(k_HH_intra)
    OVERRIDE(k_OH_intra)

    OVERRIDE(k_XH_coul)
    OVERRIDE(k_XO_coul)

    OVERRIDE(k_xlp_main)

    OVERRIDE(d_HH_intra)
    OVERRIDE(d_OH_intra)

    OVERRIDE(d_XH_coul)
    OVERRIDE(d_XO_coul)

    OVERRIDE(d_xlp_main)

//    OVERRIDE(in_plane_gamma)
//    OVERRIDE(out_of_plane_gamma)

    OVERRIDE(r2i)
    OVERRIDE(r2f)

#   undef OVERRIDE
}

//----------------------------------------------------------------------------//

void  x2b_h2o_ion_v1_common::cart_to_vars(const double* xyz, double* v, double& s) const
{
    const double* O  = xyz;
    const double* H1 = xyz + 3;
    const double* H2 = xyz + 6;

    const double* X  = xyz + 9;

    double Xa1[3], Xa2[3];

    xpoints(m_in_plane_gamma, m_out_of_plane_gamma, O, Xa1, Xa2);

//    const double d0_intra = 0.0; //should they have different values?
//    const double d0_coul = 0.0;
//    const double d0_main = 0.0;

    using x2o::distance;

    v[0]  = var_intra(m_d_HH_intra, m_k_HH_intra, distance(H1, H2));
    v[1]  = var_intra(m_d_OH_intra, m_k_OH_intra, distance( O, H1));
    v[2]  = var_intra(m_d_OH_intra, m_k_OH_intra, distance( O, H2));

    v[3]  = var_coul(m_d_XH_coul, m_k_XH_coul, distance(X, H1));
    v[4]  = var_coul(m_d_XH_coul, m_k_XH_coul, distance(X, H2));
    v[5]  = var_coul(m_d_XO_coul, m_k_XO_coul, distance(X, O));

    v[6] = var_main(m_d_xlp_main, m_k_xlp_main, distance(Xa1, X));
    v[7] = var_main(m_d_xlp_main, m_k_xlp_main, distance(Xa2, X));

    s = f_switch(distance(O, X));

	PR(s)
	PR(v[0]);
	PR(v[1]);
	PR(v[2]);
	PR(v[3]);
	PR(v[4]);
	PR(v[5]);
	PR(v[6]);
	PR(v[7]);
} //TODO check!!

//----------------------------------------------------------------------------//

void x2b_h2o_ion_v1_common::write_cdl(std::ostream& os, unsigned DEG,
                              unsigned npoly, const double* poly) const
{
    using namespace std;

    ios::fmtflags saved_flags = os.flags();

    os << setprecision(15) << scientific
       << "netcdf x2b_h2o_ion_v1 {\n"
          "  // global attributes\n"
          "  :name = \"x2b_h2o_ion_v1<" << DEG << ">\";\n";

    x2b_h2o_ion_base::as_cdl(os);

    os << "  :k_HH_intra = " << setw(22) << m_k_HH_intra << "; // A^(-1)\n"
          "  :k_OH_intra = " << setw(22) << m_k_OH_intra << "; // A^(-1)\n"
          "  :k_XH_coul = " << setw(22) << m_k_XH_coul << "; // A^(-1)\n"
          "  :k_XO_coul = " << setw(22) << m_k_XO_coul << "; // A^(-1)\n"
          "  :k_xlp_main = " << setw(22) << m_k_xlp_main << "; // A^(-1)\n"
          "  :d_HH_intra = " << setw(22) << m_d_HH_intra << "; // A^(-1)\n"
          "  :d_OH_intra = " << setw(22) << m_d_OH_intra << "; // A^(-1)\n"
          "  :d_XH_coul = " << setw(22) << m_d_XH_coul << "; // A^(-1)\n"
          "  :d_XO_coul = " << setw(22) << m_d_XO_coul << "; // A^(-1)\n"
          "  :d_xlp_main = " << setw(22) << m_d_xlp_main << "; // A^(-1)\n"
          "  :in_plane_gamma = " << setw(22) << m_in_plane_gamma << ";\n"
          "  :out_of_plane_gamma = " << setw(22) << m_out_of_plane_gamma << ";\n"
          "  :r2i = " << setw(22) << m_r2i << "; // A\n"
          "  :r2f = " << setw(22) << m_r2f << "; // A\n"
          "  dimensions:\n"
          "  poly = " << npoly << ";\n"
          "  variables:\n"
          "    double poly(poly);\n"
          "data:\n";

    os << "poly =\n";

    {
        unsigned n(0);
        for (n = 0; n + 1 < npoly; ++n)
            os << std::setw(22) << poly[n] << ", // " << n << '\n';
        os <<  std::setw(22) << poly[n] << "; // " << n << "\n\n";
    }

    os << "}\n";
    os.flags(saved_flags);
} 

//----------------------------------------------------------------------------//

void x2b_h2o_ion_v1_common::do_load_netcdf
    (const char* fn, unsigned npoly, double* poly)
{
    assert(fn);

    int rc, ncid;
    if ((rc = nc_open(fn, NC_NOWRITE, &ncid)))
        error(rc);

    from_cdf(ncid);

#   define RETRIEVE(name) \
    if ((rc = nc_get_att_double(ncid, NC_GLOBAL, #name , &m_##name))) \
        error(rc);

    RETRIEVE(k_HH_intra)
    RETRIEVE(k_OH_intra)

    RETRIEVE(k_XH_coul)
    RETRIEVE(k_XO_coul)

    RETRIEVE(k_xlp_main)

    RETRIEVE(d_HH_intra)
    RETRIEVE(d_OH_intra)

    RETRIEVE(d_XH_coul)
    RETRIEVE(d_XO_coul)

    RETRIEVE(d_xlp_main)

    RETRIEVE(in_plane_gamma)
    RETRIEVE(out_of_plane_gamma)

    RETRIEVE(r2i)
    RETRIEVE(r2f)

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

} // namespace h2o_ion

////////////////////////////////////////////////////////////////////////////////
