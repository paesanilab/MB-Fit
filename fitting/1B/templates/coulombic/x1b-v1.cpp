#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif // HAVE_CONFIG_H

#include <cmath>
#include <cassert>

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "stuff.h"
#include "x1b-v1.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int kode)
{
    std::cerr << " ** Fatal Error in x1b_v1::load_poly_dat() ** : "
              << kode << std::endl;
    std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------//

const size_t npoly = 50;

void evpoly(const double* x, double* p)
{
    p[0] = x[2] + x[1] + x[0];
    p[1] = x[5] + x[4] + x[3];

    p[2] = x[3]*x[4] + x[3]*x[5] + x[4]*x[5];
    p[3] = x[1]*x[5] + x[0]*x[3] + x[0]*x[4] + x[1]*x[3] + x[2]*x[5] + x[2]*x[4];
    p[4] = x[1]*x[4] + x[2]*x[3] + x[0]*x[5];
    p[5] = x[1]*x[1] + x[2]*x[2] + x[0]*x[0];
    p[6] = x[0]*x[1] + x[1]*x[2] + x[0]*x[2];
    p[7] = x[4]*x[4] + x[5]*x[5] + x[3]*x[3];

    p[8] = x[1]*x[1]*x[1] + x[0]*x[0]*x[0] + x[2]*x[2]*x[2];
    p[9] = x[3]*x[4]*x[4] + x[3]*x[3]*x[5] + x[3]*x[3]*x[4] + x[4]*x[4]*x[5] + x[3]*x[5]*x[5] + x[4]*x[5]*x[5];
    p[10] = x[1]*x[2]*x[3] + x[0]*x[2]*x[5] + x[0]*x[1]*x[5] + x[0]*x[1]*x[4] + x[0]*x[2]*x[3] + x[1]*x[2]*x[4];
    p[11] = x[0]*x[0]*x[1] + x[1]*x[2]*x[2] + x[0]*x[0]*x[2] + x[0]*x[2]*x[2] + x[0]*x[1]*x[1] + x[1]*x[1]*x[2];
    p[12] = x[0]*x[3]*x[3] + x[0]*x[4]*x[4] + x[2]*x[5]*x[5] + x[1]*x[3]*x[3] + x[2]*x[4]*x[4] + x[1]*x[5]*x[5];
    p[13] = x[0]*x[0]*x[5] + x[1]*x[1]*x[4] + x[2]*x[2]*x[3];
    p[14] = x[0]*x[4]*x[5] + x[2]*x[3]*x[5] + x[2]*x[3]*x[4] + x[1]*x[3]*x[4] + x[1]*x[4]*x[5] + x[0]*x[3]*x[5];
    p[15] = x[0]*x[3]*x[4] + x[1]*x[3]*x[5] + x[2]*x[4]*x[5];
    p[16] = x[3]*x[3]*x[3] + x[4]*x[4]*x[4] + x[5]*x[5]*x[5];
    p[17] = x[2]*x[2]*x[4] + x[0]*x[0]*x[3] + x[2]*x[2]*x[5] + x[1]*x[1]*x[5] + x[0]*x[0]*x[4] + x[1]*x[1]*x[3];
    p[18] = x[0]*x[2]*x[4] + x[1]*x[2]*x[5] + x[0]*x[1]*x[3];
    p[19] = x[2]*x[3]*x[3] + x[0]*x[5]*x[5] + x[1]*x[4]*x[4];
    p[20] = x[3]*x[4]*x[5];
    p[21] = x[0]*x[1]*x[2];

    p[22] = x[2]*x[3]*x[3]*x[3] + x[1]*x[4]*x[4]*x[4] + x[0]*x[5]*x[5]*x[5];
    p[23] = x[3]*x[3]*x[3]*x[3] + x[4]*x[4]*x[4]*x[4] + x[5]*x[5]*x[5]*x[5];
    p[24] = x[0]*x[0]*x[0]*x[5] + x[1]*x[1]*x[1]*x[4] + x[2]*x[2]*x[2]*x[3];
    p[25] = x[0]*x[2]*x[3]*x[5] + x[0]*x[1]*x[4]*x[5] + x[1]*x[2]*x[3]*x[4];
    p[26] = x[0]*x[1]*x[2]*x[4] + x[0]*x[1]*x[2]*x[3] + x[0]*x[1]*x[2]*x[5];
    p[27] = x[1]*x[1]*x[1]*x[1] + x[2]*x[2]*x[2]*x[2] + x[0]*x[0]*x[0]*x[0];
    p[28] = x[3]*x[3]*x[5]*x[5] + x[3]*x[3]*x[4]*x[4] + x[4]*x[4]*x[5]*x[5];
    p[29] = x[0]*x[0]*x[1]*x[5] + x[0]*x[0]*x[2]*x[5] + x[0]*x[2]*x[2]*x[3] + x[1]*x[2]*x[2]*x[3] + x[1]*x[1]*x[2]*x[4] + x[0]*x[1]*x[1]*x[4];
    p[30] = x[0]*x[2]*x[3]*x[3] + x[0]*x[2]*x[5]*x[5] + x[1]*x[2]*x[3]*x[3] + x[0]*x[1]*x[5]*x[5] + x[0]*x[1]*x[4]*x[4] + x[1]*x[2]*x[4]*x[4];
    p[31] = x[0]*x[2]*x[4]*x[5] + x[0]*x[1]*x[3]*x[5] + x[0]*x[1]*x[3]*x[4] + x[1]*x[2]*x[4]*x[5] + x[0]*x[2]*x[3]*x[4] + x[1]*x[2]*x[3]*x[5];
    p[32] = x[1]*x[5]*x[5]*x[5] + x[1]*x[3]*x[3]*x[3] + x[0]*x[4]*x[4]*x[4] + x[2]*x[5]*x[5]*x[5] + x[0]*x[3]*x[3]*x[3] + x[2]*x[4]*x[4]*x[4];
    p[33] = x[0]*x[1]*x[1]*x[1] + x[0]*x[0]*x[0]*x[2] + x[1]*x[1]*x[1]*x[2] + x[1]*x[2]*x[2]*x[2] + x[0]*x[2]*x[2]*x[2] + x[0]*x[0]*x[0]*x[1];
    p[34] = x[2]*x[3]*x[4]*x[5] + x[0]*x[3]*x[4]*x[5] + x[1]*x[3]*x[4]*x[5];
    p[35] = x[3]*x[4]*x[4]*x[4] + x[4]*x[5]*x[5]*x[5] + x[4]*x[4]*x[4]*x[5] + x[3]*x[3]*x[3]*x[4] + x[3]*x[3]*x[3]*x[5] + x[3]*x[5]*x[5]*x[5];
    p[36] = x[1]*x[4]*x[4]*x[5] + x[2]*x[3]*x[3]*x[4] + x[2]*x[3]*x[3]*x[5] + x[0]*x[4]*x[5]*x[5] + x[1]*x[3]*x[4]*x[4] + x[0]*x[3]*x[5]*x[5];
    p[37] = x[3]*x[3]*x[4]*x[5] + x[3]*x[4]*x[4]*x[5] + x[3]*x[4]*x[5]*x[5];
    p[38] = x[0]*x[3]*x[3]*x[5] + x[2]*x[3]*x[5]*x[5] + x[0]*x[4]*x[4]*x[5] + x[1]*x[4]*x[5]*x[5] + x[1]*x[3]*x[3]*x[4] + x[2]*x[3]*x[4]*x[4];
    p[39] = x[0]*x[1]*x[1]*x[5] + x[1]*x[1]*x[2]*x[3] + x[1]*x[2]*x[2]*x[4] + x[0]*x[0]*x[2]*x[3] + x[0]*x[0]*x[1]*x[4] + x[0]*x[2]*x[2]*x[5];
    p[40] = x[0]*x[1]*x[1]*x[2] + x[0]*x[0]*x[1]*x[2] + x[0]*x[1]*x[2]*x[2];
    p[41] = x[0]*x[3]*x[3]*x[4] + x[0]*x[3]*x[4]*x[4] + x[1]*x[3]*x[5]*x[5] + x[2]*x[4]*x[4]*x[5] + x[2]*x[4]*x[5]*x[5] + x[1]*x[3]*x[3]*x[5];
    p[42] = x[0]*x[0]*x[1]*x[3] + x[1]*x[1]*x[2]*x[5] + x[0]*x[1]*x[1]*x[3] + x[1]*x[2]*x[2]*x[5] + x[0]*x[0]*x[2]*x[4] + x[0]*x[2]*x[2]*x[4];
    p[43] = x[2]*x[2]*x[5]*x[5] + x[1]*x[1]*x[5]*x[5] + x[0]*x[0]*x[4]*x[4] + x[0]*x[0]*x[3]*x[3] + x[2]*x[2]*x[4]*x[4] + x[1]*x[1]*x[3]*x[3];
    p[44] = x[2]*x[2]*x[4]*x[5] + x[0]*x[0]*x[3]*x[4] + x[1]*x[1]*x[3]*x[5];
    p[45] = x[1]*x[1]*x[1]*x[5] + x[2]*x[2]*x[2]*x[4] + x[1]*x[1]*x[1]*x[3] + x[2]*x[2]*x[2]*x[5] + x[0]*x[0]*x[0]*x[4] + x[0]*x[0]*x[0]*x[3];
    p[46] = x[2]*x[2]*x[3]*x[4] + x[2]*x[2]*x[3]*x[5] + x[0]*x[0]*x[3]*x[5] + x[0]*x[0]*x[4]*x[5] + x[1]*x[1]*x[3]*x[4] + x[1]*x[1]*x[4]*x[5];
    p[47] = x[0]*x[0]*x[2]*x[2] + x[1]*x[1]*x[2]*x[2] + x[0]*x[0]*x[1]*x[1];
    p[48] = x[1]*x[1]*x[4]*x[4] + x[2]*x[2]*x[3]*x[3] + x[0]*x[0]*x[5]*x[5];
    p[49] = x[0]*x[2]*x[4]*x[4] + x[0]*x[1]*x[3]*x[3] + x[1]*x[2]*x[5]*x[5];

}

//----------------------------------------------------------------------------//

double var_main(const double& k, const double& d, const double& d0)
{
    return std::exp(-k*(d - d0))/d;
}

//----------------------------------------------------------------------------//

double var_intra(const double& k, const double& d, const double& d0)
{
    return std::exp(-k*(d - d0))/d;
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h3o {

//----------------------------------------------------------------------------//

x1b_v1::x1b_v1()
{
    m_kOH_intra = 1.0;
    m_kHH_intra = 1.0;

    m_dOH_intra = 2.0;
    m_dHH_intra = 2.0;

    m_poly = new double[npoly];
}

//----------------------------------------------------------------------------//

x1b_v1::~x1b_v1()
{
    delete[] m_poly;
}

//----------------------------------------------------------------------------//

std::string x1b_v1::name()
{
    return "x1b_v1";
}

//----------------------------------------------------------------------------//

size_t x1b_v1::nparams()
{
    return npoly;
}

//----------------------------------------------------------------------------//

void x1b_v1::set(const double* x)
{
    std::copy(x, x + npoly, m_poly);
}

//----------------------------------------------------------------------------//

double x1b_v1::value(const double xyz[12]) const
{
    double mono[npoly];

    const double Ebase = basis(xyz, mono);

    double Epoly(0);
    for (size_t n = 0; n < npoly; ++n)
        Epoly += mono[n]*m_poly[n];

    return Epoly + Ebase;
}
//----------------------------------------------------------------------------//

double x1b_v1::basis(const double xyz[12], double* mmm) const
{
    const double* Oa  = xyz;
    const double* Ha1 = xyz + 3;
    const double* Ha2 = xyz + 6;
    const double* Ha3 = xyz + 9;

    double v[6];

    v[0] = var_intra(m_kHH_intra, distance(Ha1, Ha2), m_dHH_intra);
    v[1] = var_intra(m_kHH_intra, distance(Ha1, Ha3), m_dHH_intra);
    v[2] = var_intra(m_kHH_intra, distance(Ha2, Ha3), m_dHH_intra);

    v[3] = var_intra(m_kOH_intra, distance( Oa, Ha1), m_dOH_intra);
    v[4] = var_intra(m_kOH_intra, distance( Oa, Ha2), m_dOH_intra);
    v[5] = var_intra(m_kOH_intra, distance( Oa, Ha3), m_dOH_intra);

    evpoly(v, mmm);
    return 0;
}

void x1b_v1::as_poly_dat(std::ostream& os) const
{
    using namespace std;

    std::cout << " npoly = " << npoly ;

    ios::fmtflags saved_flags = os.flags();

    os << setprecision(15) << scientific;

//    x1b_v1::as_cdl(os);

    os << "kHH_intra " << setw(22) << m_kHH_intra << " // A^(-1)\n";
    os << "kOH_intra " << setw(22) << m_kOH_intra << " // A^(-1)\n";
    os << "dHH_intra " << setw(22) << m_dHH_intra << " // A\n";
    os << "dOH_intra " << setw(22) << m_dOH_intra << " // A\n";
    os << "polynomials: \n";

    std::cout << "npoly = " << npoly << " // \n";
    unsigned n(0);
    for (n = 0; n + 1 < npoly; ++n){
        os << std::setw(22) << m_poly[n] << " " << n << " // " << '\n';
    }
    os <<  std::setw(22) << m_poly[n] << " " << n << " // " << "\n\n";

    os << "}\n";
    os.flags(saved_flags);
}

//----------------------------------------------------------------------------//

void x1b_v1::load_poly_dat(const char* fn)
{
    assert(fn);
    std::ifstream infile(fn);

    std::string line;

    bool polynomial_flag = false;

    while (std::getline(infile, line))
    {
        std::string identifier;
        std::istringstream iss(line);
        iss >> identifier;

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
	    identifier_to_poly >> m_poly[poly_num];

	    if (poly_num == npoly - 1) polynomial_flag = false;

        }

	if (polynomial_flag == false && identifier == "polynomials:")
	    polynomial_flag = true;

    } 

}

//----------------------------------------------------------------------------//

void x1b_v1::set_nonlinear_parameters(const double* xxx)
{
    m_kOH_intra = *xxx++;
    m_kHH_intra = *xxx++;

    m_dOH_intra = *xxx++;
    m_dHH_intra = *xxx++;
}

//----------------------------------------------------------------------------//

void x1b_v1::get_nonlinear_parameters(double* xxx) const
{
    *xxx++ = m_kOH_intra;
    *xxx++ = m_kHH_intra;

    *xxx++ = m_dOH_intra;
    *xxx++ = m_dHH_intra;
}

//----------------------------------------------------------------------------//

bool x1b_v1::nonlinear_parameters_out_of_range() const
{
    const double k_min = -2.0;
    const double k_max = 2.0;

    return m_dOH_intra < 0 || m_dHH_intra < 0 || m_kOH_intra < k_min ||
	   m_kOH_intra > k_max || m_kHH_intra < k_min || m_kHH_intra > k_max;
}

//----------------------------------------------------------------------------//

} // namespace x2o

////////////////////////////////////////////////////////////////////////////////
