
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <netcdf.h>

#include "stuff.h"
#include "constants.h"
#include "x2b_A1B2_A1B2_v1x.h" 
 

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int kode) {

    std::cerr << " ** Fatal Error in x2b_A1B2_A1B2_v1x::load_netcdf() ** :" 
              << nc_strerror(kode) << std::endl;
    std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------//

struct variable {
    double v_exp(const double& r0, const double& k,
                 const double* xcrd, int o, int x );

    double v_coul(const double& r0, const double& k,
      const double* xcrd, int o, int x);

    double g[3]; // diff(value, p1 - p2)
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_exp(const double& r0, const double& k,
                       const double* xcrd, int o, int x)
{
    g[0] = xcrd[o++] - xcrd[x++];
    g[1] = xcrd[o++] - xcrd[x++];
    g[2] = xcrd[o]   - xcrd[x];

    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);

    const double exp1 = std::exp(k*(r0 - r));
    const double gg = - k*exp1/r;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return exp1;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_coul(const double& r0, const double& k,
                        const double* xcrd, int o, int x)
{
    g[0] = xcrd[o++] - xcrd[x++];
    g[1] = xcrd[o++] - xcrd[x++];
    g[2] = xcrd[o]   - xcrd[x];

    const double rsq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
    const double r = std::sqrt(rsq);

    const double exp1 = std::exp(k*(r0 - r));
    const double rinv = 1.0/r;
    const double val = exp1*rinv;

    const double gg = - (k + rinv)*val*rinv;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return val;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace x2b_A1B2_A1B2 {

//----------------------------------------------------------------------------//


std::string x2b_A1B2_A1B2_v1x::name() {
    return "x2b_A1B2_A1B2_v1x";
}
void x2b_A1B2_A1B2_v1x::set_nonlinear_parameters(const double* xxx) 
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

//----------------------------------------------------------------------------//

void x2b_A1B2_A1B2_v1x::get_nonlinear_parameters(double* xxx) 
 { 
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

void x2b_A1B2_A1B2_v1x::set(const double* xxx) {
  std::copy(xxx, xxx + nparams(), m_poly);
}

//----------------------------------------------------------------------------//


void x2b_A1B2_A1B2_v1x::write_cdl(std::ostream& os, unsigned DEG,
                              unsigned npoly, const double* poly) const
{
    using namespace std;

    ios::fmtflags saved_flags = os.flags();

    os << setprecision(15) << scientific
       << "netcdf x2b_A1B2_A1B2_v1x {" << endl
       << "  // global attributes " << endl
       << "  :name = \"x2b_A1B2_A1B2_v1x<" << DEG << ">\";" << endl;

    // x2b_A1B2_A1B2_v1x::as_cdl(os);
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

         << "  :r2i = " << setw(22) << m_r2i << "; // A" << endl
         << "  :r2f = " << setw(22) << m_r2f << "; // A" << endl
         << "  dimensions:" << endl
         << "  poly = " << npoly << ";" << endl
         << "  variables:" << endl
         << "    double poly(poly);" << endl
         << "data:" << endl;

    os << "poly =" << endl;

    {
        unsigned n(0);
        for (n = 0; n + 1 < npoly; ++n)
            os << std::setw(22) << poly[n] << ", // " << n << endl;
        os <<  std::setw(22) << poly[n] << "; // " << n << endl << endl;
    }

    os << "}" << endl;
    os.flags(saved_flags);
}

//----------------------------------------------------------------------------//
bool x2b_A1B2_A1B2_v1x::nonlinear_parameters_out_of_range() const { 

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




// ##DEFINE HERE## need to think how to generalize this function...

void  x2b_A1B2_A1B2_v1x::cart_to_vars(const double* xyz, double* v, double& s) const { 

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

    double g = 0;
    s = f_switch(distance(A, A_2), g);
#define PR(x)
  PR(s);
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

//----------------------------------------------------------------------------//

void x2b_A1B2_A1B2_v1x::load_netcdf(const char* fn)
{
    assert(fn);

    int rc, ncid;
    if ((rc = nc_open(fn, NC_NOWRITE, &ncid)))
        error(rc);

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


    RETRIEVE(r2i)
    RETRIEVE(r2f)

#   undef RETRIEVE

    int varid;

    if ((rc = nc_inq_varid(ncid, "poly", &varid)))
        error(rc);

    for (size_t n = 0; n < poly_type::size; ++n) {
        if ((rc = nc_get_var1_double(ncid, varid, &n, m_poly + n)))
            error(rc);
    }

    if ((rc = nc_close(ncid)))
        error(rc);
}

//----------------------------------------------------------------------------//

double x2b_A1B2_A1B2_v1x::f_switch(const double& r, double& g) const
{
    if (r > m_r2f) {
        g = 0.0;
        return 0.0;
    } else if (r > m_r2i) {
        const double t1 = M_PI/(m_r2f - m_r2i);
        const double x = (r - m_r2i)*t1;
        g = - std::sin(x)*t1/2.0;
        return (1.0 + std::cos(x))/2.0;
    } else {
        g = 0.0;
        return 1.0;
    }
}

//----------------------------------------------------------------------------//

double x2b_A1B2_A1B2_v1x::eval(const double* mon1, const double* mon2 ) const
{
    // the switch

    // ##DEFINE HERE## right now it assumes 1st atom of each monomer
    const double d12[3] = {mon1[0] -  mon2[0],
                           mon1[1] -  mon2[1],
                           mon1[2] -  mon2[2]};

    const double r12sq = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
    const double r12 = std::sqrt(r12sq);

    if (r12 > m_r2f)
        return 0.0;

    // offsets
    const int A  = 0;
    const int B1 = 3;
    const int B2 = 6;

    const int A_2  = 9;
    const int B1_2 = 12;
    const int B2_2 = 15;


    // ##DEFINE HERE## Careful... double check before using
    double xcrd[18]; // coordinates including extra-points

    std::copy(mon1, mon1 + 9, xcrd);
    std::copy(mon2, mon2 + 9, xcrd + 9);

    // ##DEFINE HERE## For now, assuming no extra points are present (lone pairs, etc)
    // the extra-points

    // ##DEFINE HERE## Need to find a way to generalize that
    
    double v[15]; 
    
    variable ctxt[15];
std::cerr << "B1_before = " << B1 << std::endl;
    v[0]  = ctxt[0].v_exp(m_d_intra_BB, m_k_intra_BB, xcrd, B1, B2);
std::cerr << "B1_after = " << B1 << std::endl;
    v[1]  = ctxt[1].v_exp(m_d_intra_BB, m_k_intra_BB, xcrd, B1_2, B2_2);
    v[2]  = ctxt[2].v_exp(m_d_intra_AB, m_k_intra_AB, xcrd, A, B1);
    v[3]  = ctxt[3].v_exp(m_d_intra_AB, m_k_intra_AB, xcrd, A, B2);
    v[4]  = ctxt[4].v_exp(m_d_intra_AB, m_k_intra_AB, xcrd, A_2, B1_2);
    v[5]  = ctxt[5].v_exp(m_d_intra_AB, m_k_intra_AB, xcrd, A_2, B2_2);
    
    v[6]  = ctxt[6].v_exp(m_d_BB, m_k_BB, xcrd, B1, B1_2);
    v[7]  = ctxt[7].v_exp(m_d_BB, m_k_BB, xcrd, B1, B2_2);
    v[8]  = ctxt[8].v_exp(m_d_BB, m_k_BB, xcrd, B2, B1_2);
    v[9]  = ctxt[9].v_exp(m_d_BB, m_k_BB, xcrd, B2, B2_2);
    
    v[10]  = ctxt[10].v_exp(m_d_AB, m_k_AB, xcrd, A, B1_2);
    v[11]  = ctxt[11].v_exp(m_d_AB, m_k_AB, xcrd, A, B2_2);
    v[12]  = ctxt[12].v_exp(m_d_AB, m_k_AB, xcrd, A_2, B1);
    v[13]  = ctxt[13].v_exp(m_d_AB, m_k_AB, xcrd, A_2, B2);
    
    v[14]  = ctxt[14].v_exp(m_d_AA, m_k_AA, xcrd, A, A_2);
    
    const double E_poly = mb_system::poly_model::eval(m_poly, v);
    
    // the switch

    double gsw;
    const double sw = f_switch(r12, gsw);

    return sw*E_poly;
}

double x2b_A1B2_A1B2_v1x::operator()(const double crd[18]) const
{
    const double E_poly = eval(crd, crd + 9);

    return E_poly;
}

} // namespace x2b_A1B2_A1B2

////////////////////////////////////////////////////////////////////////////////
