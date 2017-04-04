#ifdef HAVE_CONFIG_H
#   include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

#include "x1b-v1x.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int errnum)
{
    std::cerr << " ** Fatal Error in x1b_v1x::load_poly_dat() ** : "
              << std::endl;
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

    return std::exp(-k*(d - r0))/d;
}

//----------------------------------------------------------------------------//

void g_intra(const double& g, const double& k, const double& r0,
             const double* a1, const double* a2,
             double* g1, double* g2)
{
    const double dx[3] = {a1[0] - a2[0],
                          a1[1] - a2[1],
                          a1[2] - a2[2]};
    const double dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const double d = std::sqrt(dsq);

    const double gg = - k*g*std::exp(-k*(d - r0))/dsq - (g*std::exp(-k*(d-r0))/(dsq*d));

    for (int i = 0; i < 3; ++i) {
        g1[i] += gg*dx[i];
        g2[i] -= gg*dx[i];
    }
}

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
    return std::exp(-k*(d - r0))/d;
}

//----------------------------------------------------------------------------//

void g_main(const double& g, const double& k, const double& r0,
            const double* a1, const double* a2,
            double* g1, double* g2)
{
    const double dx[3] = {a1[0] - a2[0],
                          a1[1] - a2[1],
                          a1[2] - a2[2]};
    const double dsq = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
    const double d = std::sqrt(dsq);

//    const double gg = g*(1.0/d - k)*std::exp(-k*(d - r0));
    const double gg = - k*g*std::exp(-k*(d - r0))/d - (g*std::exp(-k*(d-r0))/(dsq*d));

    for (int i = 0; i < 3; ++i) {
        g1[i] += gg*dx[i];
        g2[i] -= gg*dx[i];
    }
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h3o {

//----------------------------------------------------------------------------//

double x1b_v1x::operator()(const double* w1) const
{
    const double* Oa  = w1;
    const double* Ha1 = w1 + 3;
    const double* Ha2 = w1 + 6;
    const double* Ha3 = w1 + 9;

    double x[36];

    x[0] = v_intra(m_kHH_intra, m_dHH_intra, Ha1, Ha2);
    x[1] = v_intra(m_kHH_intra, m_dHH_intra, Ha1, Ha3);
    x[2] = v_intra(m_kHH_intra, m_dHH_intra, Ha2, Ha3);
    x[3] = v_intra(m_kOH_intra, m_dOH_intra,  Oa, Ha1);
    x[4] = v_intra(m_kOH_intra, m_dOH_intra,  Oa, Ha2);
    x[5] = v_intra(m_kOH_intra, m_dOH_intra,  Oa, Ha3);

    double retval = poly_1b_v1x::eval(m_coeffs, x);

    return retval;
}

//----------------------------------------------------------------------------//

double x1b_v1x::operator()(const double* w1, double* g1) const
{
    const double* Oa  = w1;
    const double* Ha1 = w1 + 3;
    const double* Ha2 = w1 + 6;
    const double* Ha3 = w1 + 9;

    double x[6];

    x[0] = v_intra(m_kHH_intra, m_dHH_intra, Ha1, Ha2);
    x[1] = v_intra(m_kHH_intra, m_dHH_intra, Ha1, Ha3);
    x[2] = v_intra(m_kHH_intra, m_dHH_intra, Ha2, Ha3);
    x[3] = v_intra(m_kOH_intra, m_dOH_intra,  Oa, Ha1);
    x[4] = v_intra(m_kOH_intra, m_dOH_intra,  Oa, Ha2);
    x[5] = v_intra(m_kOH_intra, m_dOH_intra,  Oa, Ha3);

    double g[6];
    double retval = poly_1b_v1x::eval(m_coeffs, x, g);

    double* gOa  = g1;
    double* gHa1 = g1 + 3;
    double* gHa2 = g1 + 6;
    double* gHa3 = g1 + 9;

    g_intra(g[0], m_kHH_intra, m_dHH_intra, Ha1, Ha2, gHa1, gHa2);
    g_intra(g[1], m_kHH_intra, m_dHH_intra, Ha1, Ha3, gHa1, gHa3);
    g_intra(g[2], m_kHH_intra, m_dHH_intra, Ha2, Ha3, gHa2, gHa3);
    g_intra(g[3], m_kOH_intra, m_dOH_intra,  Oa, Ha1,  gOa, gHa1);
    g_intra(g[4], m_kOH_intra, m_dOH_intra,  Oa, Ha2,  gOa, gHa2);
    g_intra(g[5], m_kOH_intra, m_dOH_intra,  Oa, Ha3,  gOa, gHa3);

    return retval;
}

//----------------------------------------------------------------------------//

void x1b_v1x::load_poly_dat(const char* fn)
{
    assert(fn);
    std::ifstream infile(fn);

    std::string line;

    bool polynomial_flag = false;

    while (std::getline(infile, line))
    {
        std::string identifier;
        std::istringstream iss(line);
        iss >>  identifier;

        if (identifier == "kHH_intra")
            iss >> m_kHH_intra;

        if (identifier == "dHH_intra")
            iss >> m_dHH_intra;

        if (identifier == "kOH_intra")
            iss >> m_kOH_intra;

        if (identifier == "dOH_intra")
            iss >> m_dOH_intra;

        if (polynomial_flag == true){
            int poly_num;
            iss >> poly_num;

            std::istringstream identifier_to_poly(identifier);
            identifier_to_poly >> m_coeffs[poly_num];
        }

        if (polynomial_flag == false && identifier == "polynomials:")
            polynomial_flag = true;

    }

}
//----------------------------------------------------------------------------//

} // namespace h3o

////////////////////////////////////////////////////////////////////////////////
