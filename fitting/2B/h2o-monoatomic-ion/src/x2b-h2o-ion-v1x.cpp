#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>

#include <netcdf.h>

#include "stuff.h"
#include "x2b-h2o-ion-v1x.h"
#include "constants.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int kode)
{
    std::cerr << " ** Fatal Error in x2b_h2o_ion_v1x::load_netcdf() ** : "
              << nc_strerror(kode) << std::endl;
    std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------//

struct variable {
    double v_exp(const double& r0, const double& k,
                 const double* xcrd, int o, int x );

    double v_coul(const double& r0, const double& k,
		  const double* xcrd, int o, int x);

//    void grads(const double& gg, double* xgrd, int o, int x) const;

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

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

/*void variable::grads(const double& gg, double* xgrd, int o, int x ) const
{
    for (int i = 0; i < 3; ++i) {
        const double d = gg*g[i];

        xgrd[o++] += d;
        xgrd[ x++] -= d;
    }
}
*/
//----------------------------------------------------------------------------//

struct monomer {
    double oh1[3];
    double oh2[3];

    void setup(const double* ohh,
               const double& in_plane_g, const double& out_of_plane_g,
               double x1[3], double x2[3]);

/*    void grads(const double* g1, const double* g2,
               const double& in_plane_g, const double& out_of_plane_g,
  	       double* grd) const;
*/
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void monomer::setup(const double* ohh,
                    const double& in_plane_g, const double& out_of_plane_g,
                    double* x1, double* x2)
{
    for (int i = 0; i < 3; ++i) {
        oh1[i] = ohh[i + 3] - ohh[i];
        oh2[i] = ohh[i + 6] - ohh[i];
    }

    const double v[3] = {
        oh1[1]*oh2[2] - oh1[2]*oh2[1],
        oh1[2]*oh2[0] - oh1[0]*oh2[2],
        oh1[0]*oh2[1] - oh1[1]*oh2[0]
    };

    for (int i = 0; i < 3; ++i) {
        const double in_plane = ohh[i] + 0.5*in_plane_g*(oh1[i] + oh2[i]);
        const double out_of_plane = out_of_plane_g*v[i];

        x1[i] = in_plane + out_of_plane;
        x2[i] = in_plane - out_of_plane;
    }
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

//don't need gradients for now TODO
/*void monomer::grads(const double* g1, const double* g2,
                    const double& in_plane_g, const double& out_of_plane_g,
                    double* grd) const
{
    const double gm[3] = {
        g1[0] - g2[0],
        g1[1] - g2[1],
        g1[2] - g2[2]
    };

    const double t1[3] = {
        oh2[1]*gm[2] - oh2[2]*gm[1],
        oh2[2]*gm[0] - oh2[0]*gm[2],
        oh2[0]*gm[1] - oh2[1]*gm[0]
    };

    const double t2[3] = {
        oh1[1]*gm[2] - oh1[2]*gm[1],
        oh1[2]*gm[0] - oh1[0]*gm[2],
        oh1[0]*gm[1] - oh1[1]*gm[0]
    };

    for (int i = 0; i < 3; ++i) {
        const double gsum = g1[i] + g2[i];
        const double in_plane = 0.5*in_plane_g*gsum;

        const double gh1 = in_plane + out_of_plane_g*t1[i];
        const double gh2 = in_plane - out_of_plane_g*t2[i];

        grd[i + 0] += gsum - (gh1 + gh2); // O
        grd[i + 3] += gh1; // H1
        grd[i + 6] += gh2; // H2
    }
}
*/
//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h2o_ion {

//----------------------------------------------------------------------------//

std::string x2b_h2o_ion_v1x::name()
{
    return "x2b_h2o_ion_v1x";
}

//----------------------------------------------------------------------------//

void x2b_h2o_ion_v1x::load_netcdf(const char* fn)
{
    assert(fn);

    int rc, ncid;
    if ((rc = nc_open(fn, NC_NOWRITE, &ncid)))
        error(rc);

    x2b_h2o_ion_base::from_cdf(ncid);

#   define RETRIEVE(name) \
    if ((rc = nc_get_att_double(ncid, NC_GLOBAL, #name , &m_##name))) \
        error(rc);

    RETRIEVE(k_HH_intra)
    RETRIEVE(k_OH_intra)

    RETRIEVE(k_XH_coul)
    RETRIEVE(k_XO_coul)

    RETRIEVE(k_xlpmain)

    RETRIEVE(in_plane_gamma)
    RETRIEVE(out_of_plane_gamma)

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

double x2b_h2o_ion_v1x::f_switch(const double& r, double& g) const
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

double x2b_h2o_ion_v1x::polynomial(const double* w1, const double* x ) const
{
    // the switch

    const double dOX[3] = {w1[0] -  x[0],
                           w1[1] -  x[1],
                           w1[2] -  x[2]};

    const double rOXsq = dOX[0]*dOX[0] + dOX[1]*dOX[1] + dOX[2]*dOX[2];
    const double rOX = std::sqrt(rOXsq);

    if (rOX > m_r2f)
        return 0.0;

    // offsets

    const int O  = 0;
    const int H1 = 3;
    const int H2 = 6;

    const int X  = 9;
    
    const int lp1 = 12;
    const int lp2 = 15;

    double xcrd[18]; // coordinates including extra-points

    std::copy(w1, w1 + 9, xcrd);
    std::copy( x,  x + 3, xcrd + 9);

    // the extra-points

    monomer m; 

    m.setup(xcrd + O,
             m_in_plane_gamma, m_out_of_plane_gamma,
             xcrd + lp1, xcrd + lp2);

//    mb.setup(xcrd + Ob,
//             m_in_plane_gamma, m_out_of_plane_gamma,
//             xcrd + Xb1, xcrd + Xb2);

    // variables

    //whats up with these?
    const double d0_intra = 1.0;
    const double d0_inter = 4.0; 

    double v[8]; // stored separately (gets passed to poly::eval)

    variable ctxt[8];

    v[0] = ctxt[0].v_exp(d0_intra, m_k_HH_intra, xcrd, H1, H2);

    v[1] = ctxt[1].v_exp(d0_intra, m_k_OH_intra, xcrd, O, H1);
    v[2] = ctxt[2].v_exp(d0_intra, m_k_OH_intra, xcrd, O, H2);

    v[3] = ctxt[3].v_coul(d0_inter, m_k_XH_coul, xcrd, X, H1);
    v[4] = ctxt[4].v_coul(d0_inter, m_k_XH_coul, xcrd, X, H2);

    v[5] = ctxt[5].v_coul(d0_inter, m_k_XO_coul, xcrd, X, O);

    v[6] = ctxt[6].v_exp(d0_inter, m_k_xlp_main, xcrd, X, lp1);
    v[7] = ctxt[7].v_exp(d0_inter, m_k_xlp_main, xcrd, X, lp2);

    const double E_poly = poly_2b_h2o_ion_v1x::eval(m_poly, v);

    // the switch

    double gsw;
    const double sw = f_switch(rOX, gsw);

    return sw*E_poly;
}

//----------------------------------------------------------------------------//

/*double x2b_h2o_ion_v1x::polynomial
    (const double* w1, const double*  x, double* g1, double* g2) const
{
    // the switch

    const double dOX[3] = {w1[0] - x[0],
                           w1[1] - x[1],
                           w1[2] - x[2]};

    const double rOXsq = dOX[0]*dOX[0] + dOX[1]*dOX[1] + dOX[2]*dOX[2];
    const double rOX = std::sqrt(rOXsq);

    if (rOX > m_r2f)
        return 0.0;

    // offsets

    const int O  = 0;
    const int H1 = 3;
    const int H2 = 6;

    const int X  = 9;

    const int lp1 = 12;
    const int lp2 = 15;

    double xcrd[18]; // coordinates including extra-points

    std::copy(w1, w1 + 9, xcrd);
    std::copy(x ,  x + 3, xcrd + 9);

    // the extra-points

    monomer m;

    m.setup(xcrd + O,
             m_in_plane_gamma, m_out_of_plane_gamma,
             xcrd + lp1, xcrd + lp2);

//    mb.setup(xcrd + Ob,
//             m_in_plane_gamma, m_out_of_plane_gamma,
//             xcrd + Xb1, xcrd + Xb2);

    // variables

    const double d0_intra = 1.0;
    const double d0_inter = 4.0;

    double v[8]; // stored separately (gets passed to poly::eval)

    variable ctxt[8];

    v[0] = ctxt[0].v_exp(d0_intra, m_k_HH_intra, xcrd, H1, H2);

    v[1] = ctxt[1].v_exp(d0_intra, m_k_OH_intra, xcrd, O, H1);
    v[2] = ctxt[2].v_exp(d0_intra, m_k_OH_intra, xcrd, O, H2);

    v[3] = ctxt[3].v_coul(d0_inter, m_k_XH_coul, xcrd, X, H1);
    v[4] = ctxt[4].v_coul(d0_inter, m_k_XH_coul, xcrd, X, H2);

    v[5] = ctxt[5].v_coul(d0_inter, m_k_XO_coul, xcrd, X, O);

    v[6] = ctxt[6].v_exp(d0_inter, m_k_xlp_main, xcrd, X, lp1);
    v[7] = ctxt[7].v_exp(d0_inter, m_k_xlp_main, xcrd, X, lp2);

    double g[19];
    const double E_poly = poly_2b_h2o_ion_v1x::eval(m_poly, v, g);

    double xgrd[18];
    std::fill(xgrd, xgrd + 18, 0.0);

    ctxt[0].grads(g[0], xgrd, H1, H2);

    ctxt[1].grads(g[1], xgrd, O, H1);
    ctxt[2].grads(g[2], xgrd, O, H2);

    ctxt[3].grads(g[3], xgrd, X, H1);
    ctxt[4].grads(g[4], xgrd, X, H2);

    ctxt[5].grads(g[5], xgrd, X, O);

    ctxt[6].grads(g[7], xgrd, X, lp1);
    ctxt[7].grads(g[7], xgrd, X, lp2);

    // distribute gradients w.r.t. the X-points

    m.grads(xgrd + lp1, xgrd + lp2,
             m_in_plane_gamma, m_out_of_plane_gamma,
             xgrd + O);

//    mb.grads(xgrd + Xb1, xgrd + Xb2,
//             m_in_plane_gamma, m_out_of_plane_gamma,
//             xgrd + Ob);

    // the switch

    double gsw;
    const double sw = f_switch(rOX, gsw);

    for (int i = 0; i < 9; ++i) {
        g1[i] += sw*xgrd[i];
//        g2[i] += sw*xgrd[i + 9]; what happens to the chloride? FIXME
    }

    // gradient of the switch

    gsw *= E_poly/rOX;
    for (int i = 0; i < 3; ++i) {
        const double d = gsw*dOX;
        g1[i] += d;
//        g2[i] -= d;
    }

    return sw*E_poly;
}
*/
//----------------------------------------------------------------------------//
//TODO removed all the gradient parts modified this for no gradients 
double x2b_h2o_ion_v1x::operator()(const double crd[18], double grd[18]) const
{
    const double E_base = x2b_base::value(crd);
    const double E_poly = polynomial(crd, crd + 9);

    return E_base + E_poly;
}

//----------------------------------------------------------------------------//

} // namespace h2o_ion

////////////////////////////////////////////////////////////////////////////////
