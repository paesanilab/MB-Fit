#ifdef HAVE_CONFIG_H
#   include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <iostream>

#include <netcdf.h>

#include "x3b-h2o-ion-v1x.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int errnum)
{
    std::cerr << " ** Fatal Error in x3b_h2o_ion_v1x::load_netcdf() ** : "
              << nc_strerror(errnum) << std::endl;
    std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------//

double v_intra(const double& k, const double& r0,
               const double* a1, const double* a2)
{
    const double dx[3] = {a1[0] - a2[0],
                          a1[1] - a2[1],
                          a1[2] - a2[2]};
    const double dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const double d = std::sqrt(dsq);

    return std::exp(-k*(d - r0));
}

//----------------------------------------------------------------------------//

/*void g_intra(const double& g, const double& k, const double& r0,
             const double* a1, const double* a2,
             double* g1, double* g2)
{
    const double dx[3] = {a1[0] - a2[0],
                          a1[1] - a2[1],
                          a1[2] - a2[2]};
    const double dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const double d = std::sqrt(dsq);

    const double gg = - k*g*std::exp(-k*(d - r0))/d;

    for (int i = 0; i < 3; ++i) {
        g1[i] += gg*dx[i];
        g2[i] -= gg*dx[i];
    }
}
*/
//----------------------------------------------------------------------------//

double v_main(const double& k, const double& r0,
              const double* a1, const double* a2)
{
    const double dx[3] = {a1[0] - a2[0],
                          a1[1] - a2[1],
                          a1[2] - a2[2]};
    const double dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const double d = std::sqrt(dsq);

//    return std::exp(-k*(d - r0))*d;
    return std::exp(-k*(d - r0));
}

//----------------------------------------------------------------------------//

/*void g_main(const double& g, const double& k, const double& r0,
            const double* a1, const double* a2,
            double* g1, double* g2)
{
    const double dx[3] = {a1[0] - a2[0],
                          a1[1] - a2[1],
                          a1[2] - a2[2]};
    const double dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const double d = std::sqrt(dsq);

//    const double gg = g*(1.0/d - k)*std::exp(-k*(d - r0));
    const double gg = - k*g*std::exp(-k*(d - r0))/d;

    for (int i = 0; i < 3; ++i) {
        g1[i] += gg*dx[i];
        g2[i] -= gg*dx[i];
    }
}
*/
//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h2o_ion {

//----------------------------------------------------------------------------//

double x3b_h2o_ion_v1x::f_switch(const double& r, double& g) const
{
    if (r > m_r3f) {
        g = 0.0;
        return 0.0;
    } else if (r > m_r3i) {
        const double t1 = M_PI/(m_r3f - m_r3i);
        const double x = (r - m_r3i)*t1;
        g = - std::sin(x)*t1/2.0;
        return (1.0 + std::cos(x))/2.0;
    } else {
        g = 0.0;
        return 1.0;
    }
}

//----------------------------------------------------------------------------//

double x3b_h2o_ion_v1x::operator()(const double* w1,
                           const double* w2,
                           const double* x) const
{
    const double* Oa  = w1;
    const double* Ha1 = w1 + 3;
    const double* Ha2 = w1 + 6;

    const double* Ob  = w2;
    const double* Hb1 = w2 + 3;
    const double* Hb2 = w2 + 6;

    const double* X  = x;

    double x[36];

    x[0] = v_intra(m_kHH_intra, m_dHH_intra, Ha1, Ha2);
    x[1] = v_intra(m_kHH_intra, m_dHH_intra, Hb1, Hb2);
    x[2] = v_intra(m_kOH_intra, m_dOH_intra,  Oa, Ha1);
    x[3] = v_intra(m_kOH_intra, m_dOH_intra,  Oa, Ha2);
    x[4] = v_intra(m_kOH_intra, m_dOH_intra,  Ob, Hb1);
    x[5] = v_intra(m_kOH_intra, m_dOH_intra,  Ob, Hb2);

    x[6] =  v_main(m_kHH, m_dHH, Ha1, Hb1);
    x[7] = v_main(m_kHH, m_dHH, Ha1, Hb2);
    x[8] = v_main(m_kHH, m_dHH, Ha2, Hb1);
    x[9] = v_main(m_kHH, m_dHH, Ha2, Hb2);
    x[10] = v_main(m_kOH, m_dOH,  Oa, Hb1);
    x[11] = v_main(m_kOH, m_dOH,  Oa, Hb2);
    x[12] = v_main(m_kOH, m_dOH,  Ob, Ha1);
    x[13] = v_main(m_kOH, m_dOH,  Ob, Ha2);
    x[14] = v_main(m_kOO, m_dOO,  Oa,  Ob);
    x[15] = v_main(m_kXH, m_dXH,  X,  Ha1);
    x[16] = v_main(m_kXH, m_dXH,  X,  Ha2);
    x[17] = v_main(m_kXH, m_dXH,  X,  Hb1);
    x[18] = v_main(m_kXH, m_dXH,  X,  Hb2);
    x[19] = v_main(m_kXO, m_dXO,  X,  Oa);
    x[20] = v_main(m_kXO, m_dXO,  X,  Ob);


    double retval = poly_3b_h2o_ion_v1x::eval(m_coeffs, x);

    double rab[3], rac[3], rbc[3];
    double drab(0), drac(0), drbc(0);

    for (int n = 0; n < 3; ++n) {
        rab[n] = Oa[n] - Ob[n];
        drab += rab[n]*rab[n];

        rac[n] = Oa[n] - X[n];
        drac += rac[n]*rac[n];

        rbc[n] = Ob[n] - X[n];
        drbc += rbc[n]*rbc[n];
    }

    drab = std::sqrt(drab);
    drac = std::sqrt(drac);
    drbc = std::sqrt(drbc);

    double gab, gac, gbc;

    const double sab = f_switch(drab, gab);
    const double sac = f_switch(drac, gac);
    const double sbc = f_switch(drbc, gbc);

    const double s = sab*sac + sab*sbc + sac*sbc;

    return s*retval;
}

//----------------------------------------------------------------------------//

/*double x3b_h2o_ion_v1x::operator()(const double* w1,
                           const double* w2,
                           const double* x,
                           double* g1,
                           double* g2,
                           double* g3) const
{
    const double* Oa  = w1;
    const double* Ha1 = w1 + 3;
    const double* Ha2 = w1 + 6;

    const double* Ob  = w2;
    const double* Hb1 = w2 + 3;
    const double* Hb2 = w2 + 6;

    const double* Oc  = w3;
    const double* Hc1 = w3 + 3;
    const double* Hc2 = w3 + 6;

    double x[36];

    x[0] = v_intra(m_kHH_intra, m_dHH_intra, Ha1, Ha2);
    x[1] = v_intra(m_kHH_intra, m_dHH_intra, Hb1, Hb2);
    x[2] = v_intra(m_kHH_intra, m_dHH_intra, Hc1, Hc2);
    x[3] = v_intra(m_kOH_intra, m_dOH_intra,  Oa, Ha1);
    x[4] = v_intra(m_kOH_intra, m_dOH_intra,  Oa, Ha2);
    x[5] = v_intra(m_kOH_intra, m_dOH_intra,  Ob, Hb1);
    x[6] = v_intra(m_kOH_intra, m_dOH_intra,  Ob, Hb2);
    x[7] = v_intra(m_kOH_intra, m_dOH_intra,  Oc, Hc1);
    x[8] = v_intra(m_kOH_intra, m_dOH_intra,  Oc, Hc2);

    x[9] =  v_main(m_kHH, m_dHH, Ha1, Hb1);
    x[10] = v_main(m_kHH, m_dHH, Ha1, Hb2);
    x[11] = v_main(m_kHH, m_dHH, Ha1, Hc1);
    x[12] = v_main(m_kHH, m_dHH, Ha1, Hc2);
    x[13] = v_main(m_kHH, m_dHH, Ha2, Hb1);
    x[14] = v_main(m_kHH, m_dHH, Ha2, Hb2);
    x[15] = v_main(m_kHH, m_dHH, Ha2, Hc1);
    x[16] = v_main(m_kHH, m_dHH, Ha2, Hc2);
    x[17] = v_main(m_kHH, m_dHH, Hb1, Hc1);
    x[18] = v_main(m_kHH, m_dHH, Hb1, Hc2);
    x[19] = v_main(m_kHH, m_dHH, Hb2, Hc1);
    x[20] = v_main(m_kHH, m_dHH, Hb2, Hc2);
    x[21] = v_main(m_kOH, m_dOH,  Oa, Hb1);
    x[22] = v_main(m_kOH, m_dOH,  Oa, Hb2);
    x[23] = v_main(m_kOH, m_dOH,  Oa, Hc1);
    x[24] = v_main(m_kOH, m_dOH,  Oa, Hc2);
    x[25] = v_main(m_kOH, m_dOH,  Ob, Ha1);
    x[26] = v_main(m_kOH, m_dOH,  Ob, Ha2);
    x[27] = v_main(m_kOH, m_dOH,  Ob, Hc1);
    x[28] = v_main(m_kOH, m_dOH,  Ob, Hc2);
    x[29] = v_main(m_kOH, m_dOH,  Oc, Ha1);
    x[30] = v_main(m_kOH, m_dOH,  Oc, Ha2);
    x[31] = v_main(m_kOH, m_dOH,  Oc, Hb1);
    x[32] = v_main(m_kOH, m_dOH,  Oc, Hb2);
    x[33] = v_main(m_kOO, m_dOO,  Oa,  Ob);
    x[34] = v_main(m_kOO, m_dOO,  Oa,  Oc);
    x[35] = v_main(m_kOO, m_dOO,  Ob,  Oc);

    double g[36];
    double retval = poly_3b_v2x::eval(m_coeffs, x, g);

    double rab[3], rac[3], rbc[3];
    double drab(0), drac(0), drbc(0);

    for (int n = 0; n < 3; ++n) {
        rab[n] = Oa[n] - Ob[n];
        drab += rab[n]*rab[n];

        rac[n] = Oa[n] - Oc[n];
        drac += rac[n]*rac[n];

        rbc[n] = Ob[n] - Oc[n];
        drbc += rbc[n]*rbc[n];
    }

    drab = std::sqrt(drab);
    drac = std::sqrt(drac);
    drbc = std::sqrt(drbc);

    double gab, gac, gbc;

    const double sab = f_switch(drab, gab);
    const double sac = f_switch(drac, gac);
    const double sbc = f_switch(drbc, gbc);

    const double s = sab*sac + sab*sbc + sac*sbc;

    for (int n = 0; n < 36; ++n)
        g[n] *= s;

    double* gOa  = g1;
    double* gHa1 = g1 + 3;
    double* gHa2 = g1 + 6;

    double* gOb  = g2;
    double* gHb1 = g2 + 3;
    double* gHb2 = g2 + 6;

    double* gOc  = g3;
    double* gHc1 = g3 + 3;
    double* gHc2 = g3 + 6;

    g_intra(g[0], m_kHH_intra, m_dHH_intra, Ha1, Ha2, gHa1, gHa2);
    g_intra(g[1], m_kHH_intra, m_dHH_intra, Hb1, Hb2, gHb1, gHb2);
    g_intra(g[2], m_kHH_intra, m_dHH_intra, Hc1, Hc2, gHc1, gHc2);
    g_intra(g[3], m_kOH_intra, m_dOH_intra,  Oa, Ha1,  gOa, gHa1);
    g_intra(g[4], m_kOH_intra, m_dOH_intra,  Oa, Ha2,  gOa, gHa2);
    g_intra(g[5], m_kOH_intra, m_dOH_intra,  Ob, Hb1,  gOb, gHb1);
    g_intra(g[6], m_kOH_intra, m_dOH_intra,  Ob, Hb2,  gOb, gHb2);
    g_intra(g[7], m_kOH_intra, m_dOH_intra,  Oc, Hc1,  gOc, gHc1);
    g_intra(g[8], m_kOH_intra, m_dOH_intra,  Oc, Hc2,  gOc, gHc2);

    g_main(g[9],  m_kHH, m_dHH, Ha1, Hb1, gHa1, gHb1);
    g_main(g[10], m_kHH, m_dHH, Ha1, Hb2, gHa1, gHb2);
    g_main(g[11], m_kHH, m_dHH, Ha1, Hc1, gHa1, gHc1);
    g_main(g[12], m_kHH, m_dHH, Ha1, Hc2, gHa1, gHc2);
    g_main(g[13], m_kHH, m_dHH, Ha2, Hb1, gHa2, gHb1);
    g_main(g[14], m_kHH, m_dHH, Ha2, Hb2, gHa2, gHb2);
    g_main(g[15], m_kHH, m_dHH, Ha2, Hc1, gHa2, gHc1);
    g_main(g[16], m_kHH, m_dHH, Ha2, Hc2, gHa2, gHc2);
    g_main(g[17], m_kHH, m_dHH, Hb1, Hc1, gHb1, gHc1);
    g_main(g[18], m_kHH, m_dHH, Hb1, Hc2, gHb1, gHc2);
    g_main(g[19], m_kHH, m_dHH, Hb2, Hc1, gHb2, gHc1);
    g_main(g[20], m_kHH, m_dHH, Hb2, Hc2, gHb2, gHc2);
    g_main(g[21], m_kOH, m_dOH,  Oa, Hb1,  gOa, gHb1);
    g_main(g[22], m_kOH, m_dOH,  Oa, Hb2,  gOa, gHb2);
    g_main(g[23], m_kOH, m_dOH,  Oa, Hc1,  gOa, gHc1);
    g_main(g[24], m_kOH, m_dOH,  Oa, Hc2,  gOa, gHc2);
    g_main(g[25], m_kOH, m_dOH,  Ob, Ha1,  gOb, gHa1);
    g_main(g[26], m_kOH, m_dOH,  Ob, Ha2,  gOb, gHa2);
    g_main(g[27], m_kOH, m_dOH,  Ob, Hc1,  gOb, gHc1);
    g_main(g[28], m_kOH, m_dOH,  Ob, Hc2,  gOb, gHc2);
    g_main(g[29], m_kOH, m_dOH,  Oc, Ha1,  gOc, gHa1);
    g_main(g[30], m_kOH, m_dOH,  Oc, Ha2,  gOc, gHa2);
    g_main(g[31], m_kOH, m_dOH,  Oc, Hb1,  gOc, gHb1);
    g_main(g[32], m_kOH, m_dOH,  Oc, Hb2,  gOc, gHb2);
    g_main(g[33], m_kOO, m_dOO,  Oa,  Ob,  gOa,  gOb);
    g_main(g[34], m_kOO, m_dOO,  Oa,  Oc,  gOa,  gOc);
    g_main(g[35], m_kOO, m_dOO,  Ob,  Oc,  gOb,  gOc);

    // gradients of the switching function

    gab *= (sac + sbc)*retval/drab;
    gac *= (sab + sbc)*retval/drac;
    gbc *= (sab + sac)*retval/drbc;

    retval *= s;

    for (int n = 0; n < 3; ++n) {
        gOa[n] += gab*rab[n] + gac*rac[n];
        gOb[n] += gbc*rbc[n] - gab*rab[n];
        gOc[n] -= gac*rac[n] + gbc*rbc[n];
    }

    return retval;
}
*/
//----------------------------------------------------------------------------//

void x3b_h2o_ion_v1x::load_netcdf(const char* fn)
{
    assert(fn && *fn);

    int rc, ncid;
    if ((rc = nc_open(fn, NC_NOWRITE, &ncid)))
        error(rc);

#   define RETRIEVE(name) \
    if ((rc = nc_get_att_double(ncid, NC_GLOBAL, #name , &m_##name))) \
        error(rc);

    RETRIEVE(r3i)
    RETRIEVE(r3f)

    RETRIEVE(kOH_intra)
    RETRIEVE(kHH_intra)

    RETRIEVE(kOO)
    RETRIEVE(kOH)
    RETRIEVE(kHH)

    RETRIEVE(kXH)
    RETRIEVE(kXO)
    
    RETRIEVE(dOH_intra)
    RETRIEVE(dHH_intra)

    RETRIEVE(dOO)
    RETRIEVE(dOH)
    RETRIEVE(dHH)

    RETRIEVE(dXO)
    RETRIEVE(dXH)

#   undef RETRIEVE

    int varid;

    if ((rc = nc_inq_varid(ncid, "poly3", &varid)))
        error(rc);

    for (size_t n = 0; n < ncoeffs; ++n) {
        if ((rc = nc_get_var1_double(ncid, varid, &n, m_coeffs + n)))
            error(rc);
    }

    if ((rc = nc_close(ncid)))
        error(rc);
}

//----------------------------------------------------------------------------//

} // namespace h2o_ion

////////////////////////////////////////////////////////////////////////////////
