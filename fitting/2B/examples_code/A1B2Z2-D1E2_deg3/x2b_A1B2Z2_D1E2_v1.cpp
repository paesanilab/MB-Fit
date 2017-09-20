
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iomanip>

#include <netcdf.h>

#include "stuff.h"
#include "constants.h"
#include "x2b_A1B2Z2_D1E2_v1.h" 
 

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int kode) {

    std::cerr << " ** Fatal Error in x2b_A1B2Z2_D1E2_v1::load_netcdf() ** :" 
              << nc_strerror(kode) << std::endl;
    std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------//

struct variable {
    double v_exp(const double& r0, const double& k,
                 const double * p1, const double * p2 );

    double v_coul(const double& r0, const double& k,
                  const double * p1, const double * p2 );

    double g[3]; // diff(value, p1 - p2)
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_exp(const double& r0, const double& k,
                       const double * p1, const double * p2)
{
    g[0] = p1[0] - p2[0];
    g[1] = p1[1] - p2[1];
    g[2] = p1[2] - p2[2];

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
                        const double * p1, const double * p2)
{
    g[0] = p1[0] - p2[0];
    g[1] = p1[1] - p2[1];
    g[2] = p1[2] - p2[2];

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

struct monomer {
    double oh1[3];
    double oh2[3];

    void setup(const double* ohh,
               const double& in_plane_g, const double& out_of_plane_g,
               double x1[3], double x2[3]);

    void grads(const double* g1, const double* g2,
               const double& in_plane_g, const double& out_of_plane_g,
               double* grd) const;
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

void monomer::grads(const double* g1, const double* g2,
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

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

//struct vsites {
//    //void TwoParticleAverageSite() {}
//    //void ThreeParticleAverageSite() {}
//    void OutOfPlaneSite(const double& w12, const double& w13,
//                        const double& wcross, const double x1[3],
//                        const double y1[3], const double y2[3],
//                        double vs[3]);
//    //void LocalCoordinatesSite{}
//};
//
//void vsites::OutOfPlaneSite(const double& w12,
//                            const double& w13,
//                            const double& wcross,
//                            const double x1[3],
//                            const double y1[3],
//                            const double y2[3],
//                            double vs[3]) {
//    double r12[3], r13[3];
//
//    for (int i = 0; i < 3; ++i) {
//        r12[i] = y1[i] - x1[i];
//        r13[i] = y2[i] - x1[i];
//    }
//                            
//    double rc[3];
//    rc[0] = r12[1]*r13[2] - r12[2]*r13[1];
//    rc[1] = r12[2]*r13[0] - r12[0]*r13[2];
//    rc[2] = r12[0]*r13[1] - r12[1]*r13[0];
//    
//    vs[0] = x1[0] + w12 * r12[0] + w13 * r13[0] + wcross * rc[0];
//    vs[1] = x1[1] + w12 * r12[1] + w13 * r13[1] + wcross * rc[1];
//    vs[2] = x1[2] + w12 * r12[2] + w13 * r13[2] + wcross * rc[2];
//}

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace x2b_A1B2Z2_D1E2 {

//----------------------------------------------------------------------------//


std::string x2b_A1B2Z2_D1E2_v1x::name() {
    return "x2b_A1B2Z2_D1E2_v1x";
}
void x2b_A1B2Z2_D1E2_v1x::set_nonlinear_parameters(const double* xxx) 
 { 
     m_d_intra_AB = *xxx++; 
    m_k_intra_AB = *xxx++; 
    m_d_intra_BB = *xxx++; 
    m_k_intra_BB = *xxx++; 
    m_d_intra_DE = *xxx++; 
    m_k_intra_DE = *xxx++; 
    m_d_intra_EE = *xxx++; 
    m_k_intra_EE = *xxx++; 
    m_d_AD = *xxx++; 
    m_k_AD = *xxx++; 
    m_d_AE = *xxx++; 
    m_k_AE = *xxx++; 
    m_d_BD = *xxx++; 
    m_k_BD = *xxx++; 
    m_d_BE = *xxx++; 
    m_k_BE = *xxx++; 
    m_d_DZ = *xxx++; 
    m_k_DZ = *xxx++; 
    m_d_EZ = *xxx++; 
    m_k_EZ = *xxx++; 

}

//----------------------------------------------------------------------------//

void x2b_A1B2Z2_D1E2_v1x::get_nonlinear_parameters(double* xxx) 
 { 
    *xxx++ = m_d_intra_AB; 
   *xxx++ = m_k_intra_AB; 
   *xxx++ = m_d_intra_BB; 
   *xxx++ = m_k_intra_BB; 
   *xxx++ = m_d_intra_DE; 
   *xxx++ = m_k_intra_DE; 
   *xxx++ = m_d_intra_EE; 
   *xxx++ = m_k_intra_EE; 
   *xxx++ = m_d_AD; 
   *xxx++ = m_k_AD; 
   *xxx++ = m_d_AE; 
   *xxx++ = m_k_AE; 
   *xxx++ = m_d_BD; 
   *xxx++ = m_k_BD; 
   *xxx++ = m_d_BE; 
   *xxx++ = m_k_BE; 
   *xxx++ = m_d_DZ; 
   *xxx++ = m_k_DZ; 
   *xxx++ = m_d_EZ; 
   *xxx++ = m_k_EZ; 

}

//----------------------------------------------------------------------------//

void x2b_A1B2Z2_D1E2_v1x::set(const double* xxx) {
  std::copy(xxx, xxx + nparams(), m_poly);
}

//----------------------------------------------------------------------------//


void x2b_A1B2Z2_D1E2_v1x::write_cdl(std::ostream& os, unsigned DEG,
                              unsigned npoly, const double* poly) const
{
    using namespace std;

    ios::fmtflags saved_flags = os.flags();

    os << setprecision(15) << scientific
       << "netcdf x2b_A1B2Z2_D1E2_v1x {" << endl
       << "  // global attributes " << endl
       << "  :name = \"x2b_A1B2Z2_D1E2_v1x<" << DEG << ">\";" << endl;

    // x2b_A1B2Z2_D1E2_v1x::as_cdl(os);
    os        << "  :d_intra_AB = " << setw(22) << m_d_intra_AB << "; // A^(-1))" << endl 
       << "  :k_intra_AB = " << setw(22) << m_k_intra_AB << "; // A^(-1))" << endl 
       << "  :d_intra_BB = " << setw(22) << m_d_intra_BB << "; // A^(-1))" << endl 
       << "  :k_intra_BB = " << setw(22) << m_k_intra_BB << "; // A^(-1))" << endl 
       << "  :d_intra_DE = " << setw(22) << m_d_intra_DE << "; // A^(-1))" << endl 
       << "  :k_intra_DE = " << setw(22) << m_k_intra_DE << "; // A^(-1))" << endl 
       << "  :d_intra_EE = " << setw(22) << m_d_intra_EE << "; // A^(-1))" << endl 
       << "  :k_intra_EE = " << setw(22) << m_k_intra_EE << "; // A^(-1))" << endl 
       << "  :d_AD = " << setw(22) << m_d_AD << "; // A^(-1))" << endl 
       << "  :k_AD = " << setw(22) << m_k_AD << "; // A^(-1))" << endl 
       << "  :d_AE = " << setw(22) << m_d_AE << "; // A^(-1))" << endl 
       << "  :k_AE = " << setw(22) << m_k_AE << "; // A^(-1))" << endl 
       << "  :d_BD = " << setw(22) << m_d_BD << "; // A^(-1))" << endl 
       << "  :k_BD = " << setw(22) << m_k_BD << "; // A^(-1))" << endl 
       << "  :d_BE = " << setw(22) << m_d_BE << "; // A^(-1))" << endl 
       << "  :k_BE = " << setw(22) << m_k_BE << "; // A^(-1))" << endl 
       << "  :d_DZ = " << setw(22) << m_d_DZ << "; // A^(-1))" << endl 
       << "  :k_DZ = " << setw(22) << m_k_DZ << "; // A^(-1))" << endl 
       << "  :d_EZ = " << setw(22) << m_d_EZ << "; // A^(-1))" << endl 
       << "  :k_EZ = " << setw(22) << m_k_EZ << "; // A^(-1))" << endl 

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
bool x2b_A1B2Z2_D1E2_v1x::nonlinear_parameters_out_of_range() const { 
    const double k_min =  0.0 ;
    const double k_max =  3.0 ;
    const double k_min_intra =  0.0 ;
    const double k_max_intra =  2.0 ;

    const double d_min =  0.0 ;
    const double d_max =  7.0 ;
    const double d_min_intra =  0.0 ;
    const double d_max_intra =  3.0 ;

return false
       || m_d_intra_AB < d_min_intra 
       || m_d_intra_AB > d_max_intra 
       || m_k_intra_AB < k_min_intra 
       || m_k_intra_AB > k_max_intra 
       || m_d_intra_BB < d_min_intra 
       || m_d_intra_BB > d_max_intra 
       || m_k_intra_BB < k_min_intra 
       || m_k_intra_BB > k_max_intra 
       || m_d_intra_DE < d_min_intra 
       || m_d_intra_DE > d_max_intra 
       || m_k_intra_DE < k_min_intra 
       || m_k_intra_DE > k_max_intra 
       || m_d_intra_EE < d_min_intra 
       || m_d_intra_EE > d_max_intra 
       || m_k_intra_EE < k_min_intra 
       || m_k_intra_EE > k_max_intra 
       || m_d_AD < d_min 
       || m_d_AD > d_max 
       || m_k_AD < k_min 
       || m_k_AD > k_max 
       || m_d_AE < d_min 
       || m_d_AE > d_max 
       || m_k_AE < k_min 
       || m_k_AE > k_max 
       || m_d_BD < d_min 
       || m_d_BD > d_max 
       || m_k_BD < k_min 
       || m_k_BD > k_max 
       || m_d_BE < d_min 
       || m_d_BE > d_max 
       || m_k_BE < k_min 
       || m_k_BE > k_max 
       || m_d_DZ < d_min 
       || m_d_DZ > d_max 
       || m_k_DZ < k_min 
       || m_k_DZ > k_max 
       || m_d_EZ < d_min 
       || m_d_EZ > d_max 
       || m_k_EZ < k_min 
       || m_k_EZ > k_max ; 


}


void  x2b_A1B2Z2_D1E2_v1x::cart_to_vars(const double* xyz, double* v, double& s, double& gs) const { 
    // NOTE: XYZ contains ONLY the real sites. The lone pairs etc are calculated here 
    const double* A_1_a= xyz + 0;
    const double* B_1_a= xyz + 3;
    const double* B_2_a= xyz + 6;

    const double* D_1_b= xyz + 9;
    const double* E_1_b= xyz + 12;
    const double* E_2_b= xyz + 15;
    double Z_1_a[3];
    double Z_2_a[3];


//    vsites virt;
    double w12 =     -9.721486914088159e-02;  //from MBpol
    double w13 =     -9.721486914088159e-02;
    double wcross =   9.859272078406150e-02;

    monomer m;
    
    m.setup(A_1_a, w12, wcross,
             Z_1_a, Z_2_a);
                    
    variable vr[21];
    using x2o::distance;
    
    v[0]  = vr[0].v_exp(m_d_intra_AB, m_k_intra_AB, A_1_a, B_1_a);
    v[1]  = vr[1].v_exp(m_d_intra_AB, m_k_intra_AB, A_1_a, B_2_a);
    v[2]  = vr[2].v_exp(m_d_intra_BB, m_k_intra_BB, B_1_a, B_2_a);

    v[3]  = vr[3].v_exp(m_d_intra_DE, m_k_intra_DE, D_1_b, E_1_b);
    v[4]  = vr[4].v_exp(m_d_intra_DE, m_k_intra_DE, D_1_b, E_2_b);
    v[5]  = vr[5].v_exp(m_d_intra_EE, m_k_intra_EE, E_1_b, E_2_b);

    v[6]  = vr[6].v_exp(m_d_AD, m_k_AD, A_1_a, D_1_b);
    v[7]  = vr[7].v_exp(m_d_AE, m_k_AE, A_1_a, E_1_b);
    v[8]  = vr[8].v_exp(m_d_AE, m_k_AE, A_1_a, E_2_b);

    v[9]  = vr[9].v_exp(m_d_BD, m_k_BD, B_1_a, D_1_b);
    v[10]  = vr[10].v_exp(m_d_BE, m_k_BE, B_1_a, E_1_b);
    v[11]  = vr[11].v_exp(m_d_BE, m_k_BE, B_1_a, E_2_b);

    v[12]  = vr[12].v_exp(m_d_BD, m_k_BD, B_2_a, D_1_b);
    v[13]  = vr[13].v_exp(m_d_BE, m_k_BE, B_2_a, E_1_b);
    v[14]  = vr[14].v_exp(m_d_BE, m_k_BE, B_2_a, E_2_b);

    v[15]  = vr[15].v_exp(m_d_DZ, m_k_DZ, Z_1_a, D_1_b);
    v[16]  = vr[16].v_exp(m_d_EZ, m_k_EZ, Z_1_a, E_1_b);
    v[17]  = vr[17].v_exp(m_d_EZ, m_k_EZ, Z_1_a, E_2_b);

    v[18]  = vr[18].v_exp(m_d_DZ, m_k_DZ, Z_2_a, D_1_b);
    v[19]  = vr[19].v_exp(m_d_EZ, m_k_EZ, Z_2_a, E_1_b);
    v[20]  = vr[20].v_exp(m_d_EZ, m_k_EZ, Z_2_a, E_2_b);


    s = f_switch(distance(A_1_a, D_1_b), gs);
    
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
  PR(v[15]);
  PR(v[16]);
  PR(v[17]);
  PR(v[18]);
  PR(v[19]);
  PR(v[20]);

} 

//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//

void x2b_A1B2Z2_D1E2_v1x::load_netcdf(const char* fn)
{
    assert(fn);

    int rc, ncid;
    if ((rc = nc_open(fn, NC_NOWRITE, &ncid)))
        error(rc);

#   define RETRIEVE(name) \
    if ((rc = nc_get_att_double(ncid, NC_GLOBAL, #name , &m_##name))) \
        error(rc);
    RETRIEVE(d_intra_AB)
    RETRIEVE(k_intra_AB)
    RETRIEVE(d_intra_BB)
    RETRIEVE(k_intra_BB)
    RETRIEVE(d_intra_DE)
    RETRIEVE(k_intra_DE)
    RETRIEVE(d_intra_EE)
    RETRIEVE(k_intra_EE)
    RETRIEVE(d_AD)
    RETRIEVE(k_AD)
    RETRIEVE(d_AE)
    RETRIEVE(k_AE)
    RETRIEVE(d_BD)
    RETRIEVE(k_BD)
    RETRIEVE(d_BE)
    RETRIEVE(k_BE)
    RETRIEVE(d_DZ)
    RETRIEVE(k_DZ)
    RETRIEVE(d_EZ)
    RETRIEVE(k_EZ)


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

double x2b_A1B2Z2_D1E2_v1x::f_switch(const double& r, double& g) const
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

double x2b_A1B2Z2_D1E2_v1x::eval(const double* mon1, const double* mon2 ) const
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

    double xcrd[18]; // coordinates of real sites ONLY

    std::copy(mon1, mon1 + 9, xcrd);
    std::copy(mon2, mon2 + 9, xcrd + 9);
    
    double v[21]; 
    double sw = 0.0;
    double gsw = 0.0;
    cart_to_vars(xcrd, v, sw, gsw); 
    
    const double E_poly = mb_system::poly_model::eval(m_poly, v);

    return sw*E_poly;
}

double x2b_A1B2Z2_D1E2_v1x::operator()(const double crd[18]) const
{
    const double E_poly = eval(crd, crd + 9);

    return E_poly;
}

} // namespace x2b_A1B2Z2_D1E2

////////////////////////////////////////////////////////////////////////////////
