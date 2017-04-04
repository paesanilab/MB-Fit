#include <cmath>
#include <cassert>
#include <cstdlib>

#include <iomanip>

#include <netcdf.h>

#include "kvstring.h"
#include "x2b-h2o-ion-base.h"
#include "constants.h"

#include "tang-toennies.h" //dispersion damping 

//#define CHLORIDE
//#define FLUORIDE
//#define BROMIDE
//#define IODIDE
//#define DEBUG_ENERGY
#ifdef DEBUG_ENERGY
#include <iostream>
#endif 
////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//


inline double x68(const double& C6, const double& d6,
                  const double& C8, const double& d8,
                  const double* p1, const double* p2)
{
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double tt6 = x2o::tang_toennies(6, d6*r);
    const double tt8 = x2o::tang_toennies(8, d8*r);

    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;

    return - (C6*tt6 + C8*tt8*inv_rsq)*inv_r6;
}

//----------------------------------------------------------------------------//

// convenience stuff for dispersion

const double if6 = 1.0/x2o::factorial<6>();  
const double if8 = 1.0/x2o::factorial<8>();

//----------------------------------------------------------------------------//
// dispersion with gradients

inline double x68(const double& C6, const double& d6,
                  const double& C8, const double& d8,
                  const double* p1, const double* p2,
                        double* g1,       double* g2)
{
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double d6r = d6*r;
    const double tt6 = x2o::tang_toennies(6, d6r);

    const double d8r = d8*r;
    const double tt8 = x2o::tang_toennies(8, d8r);


    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;
    const double inv_r8 = inv_r6*inv_rsq;

    const double e6 = C6*tt6*inv_r6;
    const double e8 = C8*tt8*inv_r8;

    const double grd = (6*e6 + 8*e8)*inv_rsq
        - (C6*std::pow(d6, 7)*if6*std::exp(-d6r)
        +  C8*std::pow(d8, 9)*if8*std::exp(-d8r))/r;

    g1[0] += dx*grd;
    g2[0] -= dx*grd;

    g1[1] += dy*grd;
    g2[1] -= dy*grd;

    g1[2] += dz*grd;
    g2[2] -= dz*grd;

    return - (e6 + e8);
}

//----------------------------------------------------------------------------//

void error(int ec)
{
    std::cerr << " ** Fatal Error in x2b_h2o_ion_base::load() ** : "
              << nc_strerror(ec) << std::endl;
    std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------//

} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace h2o_ion {

////////////////////////////////////////////////////////////////////////////////
// conversion factors

const double au2kcalA6 = constants::Eh_kcalmol*std::pow(constants::Bohr_A, 6);
const double au2kcalA8 = constants::Eh_kcalmol*std::pow(constants::Bohr_A, 8);
#ifdef CHLORIDE
const double x2b_h2o_ion_base::C6_XO = 746.199;   //XDM Cl params 
const double x2b_h2o_ion_base::C6_XH = 306.890;

//const double x2b_h2o_ion_base::C6_XO = 79.10814*au2kcalA6;   //CamCASP 
//const double x2b_h2o_ion_base::C6_XH = 14.91344*au2kcalA6;
#endif 
#ifdef FLUORIDE
//const double x2b_h2o_ion_base::C6_XO = 721.684;   //fluoride iTTM params  
//const double x2b_h2o_ion_base::C6_XH = 512.463;

//const double x2b_h2o_ion_base::C6_XO = 598.653;   //fluoride params after 1st iter
//const double x2b_h2o_ion_base::C6_XH = 271.492;

//const double x2b_h2o_ion_base::C6_XO = 599.085;   //fluoride params after 2nd iter
//const double x2b_h2o_ion_base::C6_XH = 271.055;

const double x2b_h2o_ion_base::C6_XO = 348.864;   //XDM F params
const double x2b_h2o_ion_base::C6_XH = 128.678;
#endif 
#ifdef BROMIDE 
//const double x2b_h2o_ion_base::C6_XO = 2451.04;   //bromide  iTTM params  
//const double x2b_h2o_ion_base::C6_XH = 980.383;

//const double x2b_h2o_ion_base::C6_XO = 460.15;   //bromide  params after 1st iter
//const double x2b_h2o_ion_base::C6_XH = 532.821;

//const double x2b_h2o_ion_base::C6_XO = 854.312;   //bromide  params after 2nd iter
//const double x2b_h2o_ion_base::C6_XH = 262.463;

const double x2b_h2o_ion_base::C6_XO = 942.65;   //bromide xdm params 
const double x2b_h2o_ion_base::C6_XH = 394.168;
#endif 
#ifdef IODIDE  
//const double x2b_h2o_ion_base::C6_XO = 6320.63;   //ioidide after1stiter
//const double x2b_h2o_ion_base::C6_XH = 747.280;
 
//const double x2b_h2o_ion_base::C6_XO = 6322.18;   //ioidide iTTM params 
//const double x2b_h2o_ion_base::C6_XH = 746.732;
 
const double x2b_h2o_ion_base::C6_XO = 1294.68;   //ioidide xdm params 
const double x2b_h2o_ion_base::C6_XH = 568.156;
#endif 
const double x2b_h2o_ion_base::C8_XO = 0.0*au2kcalA8;
const double x2b_h2o_ion_base::C8_XH = 0.0*au2kcalA8;


//----------------------------------------------------------------------------//

x2b_h2o_ion_base::x2b_h2o_ion_base()
{	
    setup("");
}

//----------------------------------------------------------------------------//

void x2b_h2o_ion_base::setup(const char* str)
{
    std::map<std::string, double> kv;
    kit::kvstring_parse(str, kv);

    do_setup(kv);
}

//----------------------------------------------------------------------------//

void x2b_h2o_ion_base::do_setup(const std::map<std::string, double>& kv)
{
    // set defaults

    // damping for dispersion
  
    m_d8_XH = 2.0;
    m_d8_XO = 2.0;
   
    // i-TTM parameters 
#ifdef FLUORIDE
    //Updated after publication. Entered on 03/12/17 
    m_d6_XH = 2.69768; // range is between 1 and 10 
    m_d6_XO = 3.58619; // range is between 1 and 10 
#ifdef POL75    
    m_d6_XH = 2.71139; // range is between 1 and 10 
    m_d6_XO = 3.60872; // range is between 1 and 10 
#endif 
#ifdef POL50    
    m_d6_XH = 2.72528; // range is between 1 and 10 
    m_d6_XO = 3.63151; // range is between 1 and 10 
#endif 
#ifdef POL25    
    m_d6_XH = 2.73939; // range is between 1 and 10 
    m_d6_XO = 3.65448; // range is between 1 and 10 
#endif 
#ifdef POL0   
    m_d6_XH = 2.75372; // range is between 1 and 10 
    m_d6_XO = 3.67764; // range is between 1 and 10 
#endif 
#ifdef POL_EFF
    m_d6_XH = 2.72911; // range is between 1 and 10 
    m_d6_XO = 3.63775; // range is between 1 and 10 
#endif 
#endif 

#ifdef CHLORIDE
    //Updated after publication. Entered on 03/12/17 
    m_d6_XH = 2.78226; // range is between 1 and 10 
    m_d6_XO = 3.27542; // range is between 1 and 10 
#ifdef POL75    
    m_d6_XH = 2.79213; // range is between 1 and 10 
    m_d6_XO = 3.29506; // range is between 1 and 10 
#endif 
#ifdef POL50    
    m_d6_XH = 2.80233; // range is between 1 and 10 
    m_d6_XO = 3.31488; // range is between 1 and 10 
#endif 
#ifdef POL0   
    m_d6_XH = 2.82384; // range is between 1 and 10 
    m_d6_XO = 3.35500; // range is between 1 and 10 
#endif 
#ifdef POL25    
    m_d6_XH = 2.81289; // range is between 1 and 10 
    m_d6_XO = 3.33488; // range is between 1 and 10 
#endif 
#ifdef POL_EFF
    m_d6_XH = 2.79872; // range is between 1 and 10 
    m_d6_XO = 3.30792; // range is between 1 and 10 
#endif 
#endif 

#ifdef BROMIDE 
    //Updated after publication. Entered on 03/12/17 
    m_d6_XH = 2.79804; // range is between 1 and 10 
    m_d6_XO = 3.05825; // range is between 1 and 10 
#ifdef POL75    
    m_d6_XH = 2.81109; // range is between 1 and 10 
    m_d6_XO = 3.07325; // range is between 1 and 10 
#endif 
#ifdef POL50    
    m_d6_XH = 2.82452; // range is between 1 and 10 
    m_d6_XO = 3.08837; // range is between 1 and 10 
#endif 
#ifdef POL25    
    m_d6_XH = 2.83833; // range is between 1 and 10 
    m_d6_XO = 3.10363; // range is between 1 and 10 
#endif 
#ifdef POL0   
    m_d6_XH = 2.85250; // range is between 1 and 10 
    m_d6_XO = 3.11898; // range is between 1 and 10 
#endif 
#ifdef POL_EFF
    m_d6_XH = 2.81788; // range is between 1 and 10 
    m_d6_XO = 3.08093; // range is between 1 and 10 
#endif 
#endif 

#ifdef IODIDE 
    //Updated after publication. Entered on 10/10/16 
    m_d6_XH = 2.79911; // range is between 1 and 10 
    m_d6_XO = 2.72314; // range is between 1 and 10 
#ifdef POL75    
    m_d6_XH = 2.81094; // range is between 1 and 10 
    m_d6_XO = 2.73340; // range is between 1 and 10 
#endif 
#ifdef POL50    
    m_d6_XH = 2.82308; // range is between 1 and 10 
    m_d6_XO = 2.74381; // range is between 1 and 10 
#endif 
#ifdef POL25    
    m_d6_XH = 2.83555; // range is between 1 and 10 
    m_d6_XO = 2.75432; // range is between 1 and 10 
#endif 
#ifdef POL0   
    m_d6_XH = 2.84838; // range is between 1 and 10 
    m_d6_XO = 2.7649; // range is between 1 and 10 
#endif 
#ifdef POL_EFF
    m_d6_XH = 2.81420; // range is between 1 and 10 
    m_d6_XO = 2.73622; // range is between 1 and 10 
#endif 
#endif 

    // override from the KV
    // reads in dispersion parameter initial guesses
        
    std::map<std::string, double>::const_iterator iter;
#   define OVERRIDE(name)       \
    iter = kv.find(#name);      \
    if (iter != kv.end())       \
        m_##name = iter->second;

//    OVERRIDE(d6_XH)
//    OVERRIDE(d6_XO)
//
//    OVERRIDE(d8_XH)
//    OVERRIDE(d8_XO)

#   undef OVERRIDE

}

//----------------------------------------------------------------------------//

void x2b_h2o_ion_base::as_cdl(std::ostream& os) const
{
    
    using namespace std;

    ios::fmtflags saved_flags = os.flags();

    os << setprecision(15) << scientific
	<< "  :C6_XH = " << setw(22) << C6_XH << "; // kcal/mol * A^6\n"
	<< "  :C6_XO = " << setw(22) << C6_XO << "; // kcal/mol * A^6\n"
	<< "  :C8_XH = " << setw(22) << C8_XH << "; // kcal/mol * A^8\n"
	<< "  :C8_XO = " << setw(22) << C8_XO << "; // kcal/mol * A^8\n"
	<< "  :d6_XH = " << setw(22) << m_d6_XH << "; // A^(-1)\n"
	<< "  :d6_XO = " << setw(22) << m_d6_XO << "; // A^(-1)\n"
	<< "  :d8_XH = " << setw(22) << m_d8_XH << "; // A^(-1)\n"
	<< "  :d8_XO = " << setw(22) << m_d8_XO << "; // A^(-1)\n";

    os.flags(saved_flags);
    
}

//----------------------------------------------------------------------------//

void x2b_h2o_ion_base::from_cdf(int ncid)
{
    // the only params common to all subclasses
    // are the dispersion damping params
    
    int rc;

#   define RETRIEVE(name) \
    if ((rc = nc_get_att_double(ncid, NC_GLOBAL, #name , &m_##name))) \
        error(rc);

    RETRIEVE(d6_XH)
    RETRIEVE(d6_XO)
    
    RETRIEVE(d8_XH)
    RETRIEVE(d8_XO)
    
#   undef RETRIEVE

}

//----------------------------------------------------------------------------//

double x2b_h2o_ion_base::dispersion(const double* crd) const
{
    const double* O  = crd;
    const double* H1 = crd + 3;
    const double* H2 = crd + 6;

    const double* X  = crd + 9;

    const double XH68 = 
	x68(C6_XH, m_d6_XH, C8_XH, m_d8_XH, X, H1)
	+ x68(C6_XH, m_d6_XH, C8_XH, m_d8_XH, X, H2);

    const double XO68 = 
	x68(C6_XO, m_d6_XO, C8_XO, m_d8_XO, X, O);
    
    return XH68 + XO68;
}

//----------------------------------------------------------------------------//

double x2b_h2o_ion_base::value(const double* crd) const
{
    //The electrostatics stuff.
    double* h2o_atmcrds = new double [9];
    std::copy(crd, (crd + 9), h2o_atmcrds);

    double*  x_atmcrds = new double [3];
    std::copy((crd + 9), (crd + 12),  x_atmcrds);
    
    x2o::ttm4 h2o1(h2o_atmcrds);
     x::ion       x1( x_atmcrds);
    std::vector<molecule*> system;
    system.push_back(&h2o1);
    system.push_back(&x1 );

    int system_nsites = 0;
    double* system_sitecrds;
    double* system_charge;
    double* system_polfac;
    double* system_atmpolar;

    for (int i = 0; i < system.size(); i++) {
        system_nsites += system[i]->get_nsites();
    } //calculates total number of sites in the system

    system_sitecrds = new double [system_nsites*3];
    system_charge = new double [system_nsites];
    system_polfac = new double [system_nsites];
    system_atmpolar = new double [system_nsites]; //allocates memory to the pointers 

    int pos = 0;
    for (int i = 0; i < system.size(); i++) {

        std::copy(system[i]->get_sitecrds(), (system[i]->get_sitecrds() + system[i]->get_nsites()*3), (system_sitecrds + pos*3));

        std::copy(system[i]->get_charges(), (system[i]->get_charges() + system[i]->get_nsites()), (system_charge + pos));

	std::copy(system[i]->get_polfacs(), (system[i]->get_polfacs() + system[i]->get_nsites()), (system_polfac + pos));

        std::copy(system[i]->get_pol(), (system[i]->get_pol() + system[i]->get_nsites()), (system_atmpolar + pos));

        pos += system[i]->get_nsites();
    }

    excluded_set_type m_excluded12;
    excluded_set_type m_excluded13;
    m_excluded12.clear();
    m_excluded13.clear();

    m_excluded12.insert(std::make_pair(0, 1));
    m_excluded12.insert(std::make_pair(0, 2));
    m_excluded12.insert(std::make_pair(0, 3));
    m_excluded13.insert(std::make_pair(1, 2));
    m_excluded12.insert(std::make_pair(1, 3));
    m_excluded12.insert(std::make_pair(2, 3));

    ttm::smear_ttm4x smr; // smearing for chloride?
    smr.m_aDD_intra_12 = 0.626;
    smr.m_aDD_intra_13 = 0.055;

    //ttm::electrostatics test;

    double ener = m_electrostatics(5, system_charge, system_polfac, system_atmpolar, system_sitecrds, m_excluded12, m_excluded13, smr, 0);

    //return Eelec + Eind + dispersion(crd);

#ifdef DEBUG_ENERGY
	std::cout << "Electrostatics = " << ener << std::endl;
	std::cout << "Dispersion = " << dispersion(crd) << std::endl;
#endif

	// FIXME: free up memory
	delete[] system_sitecrds;
	delete[] system_charge;
	delete[] system_polfac;
	delete[] system_atmpolar;
	delete[] h2o_atmcrds;
	delete[] x_atmcrds;
    
//	//for Andrea
//	double disp_grad[12];
//	std::fill(disp_grad, disp_grad + 12, 0.0); 
//	std::cout << "Dispersion energy = " << dispersion(crd, crd + 9, disp_grad, disp_grad + 9) << std::endl;
//
//	for (size_t i = 0; i < 4; i++) {
//	    std::cout << "Force on atom " << i+1 << "   " <<   
//		disp_grad[3*i] << "    " << disp_grad[3*i + 1] << "    " << disp_grad[3*i + 2] << std::endl;
//	}
//	//
    return ener + dispersion(crd);
}

//----------------------------------------------------------------------------//

//TODO check if this is neccessary
/*void x2b_h2o_ion_base::components(const double* crd,
    double& Eelec, double& Eind, double& Edisp) const
{
    m_ttm4(2, crd, Eelec, 0, Eind, 0);
    //Edisp = dispersion(crd);
    
    // 0 for now
    Edisp = 0;
}
*/

//----------------------------------------------------------------------------//

double x2b_h2o_ion_base::value(const double* crd, double* grd) const
{
    double Eind, Eelec;

    double* h2o_atmcrds = new double [9];
    std::copy(crd, (crd + 9), h2o_atmcrds);

    double*  x_atmcrds = new double [3];
    std::copy((crd + 9), (crd + 12),  x_atmcrds);

    x2o::ttm4 h2o1(h2o_atmcrds);
     x::ion       x1( x_atmcrds);
    std::vector<molecule*> system;
    system.push_back(&h2o1);
    system.push_back(&x1 );

    int system_nsites = 0;
    double* system_sitecrds;
    double* system_charge;
    double* system_polfac;
    double* system_atmpolar;

    for (int i = 0; i < system.size(); i++) {
        system_nsites += system[i]->get_nsites();
    } //calculates total number of sites in the system

    system_sitecrds = new double [system_nsites*3];
    system_charge = new double [system_nsites];
    system_polfac = new double [system_nsites];
    system_atmpolar = new double [system_nsites]; //allocates memory to the pointers 

    int pos = 0;
    for (int i = 0; i < system.size(); i++) {

        std::copy(system[i]->get_sitecrds(), (system[i]->get_sitecrds() + system[i]->get_nsites()*3), (system_sitecrds + pos*3));

        std::copy(system[i]->get_charges(), (system[i]->get_charges() + system[i]->get_nsites()), (system_charge + pos));

        std::copy(system[i]->get_polfacs(), (system[i]->get_polfacs() + system[i]->get_nsites()), (system_polfac + pos));

        std::copy(system[i]->get_pol(), (system[i]->get_pol() + system[i]->get_nsites()), (system_atmpolar + pos));

        pos += system[i]->get_nsites();
    }

    excluded_set_type m_excluded12;
    excluded_set_type m_excluded13;
    m_excluded12.clear();
    m_excluded13.clear();

    m_excluded12.insert(std::make_pair(0, 1));
    m_excluded12.insert(std::make_pair(0, 2));
    m_excluded12.insert(std::make_pair(0, 3));
    m_excluded13.insert(std::make_pair(1, 2));
    m_excluded12.insert(std::make_pair(1, 3));
    m_excluded12.insert(std::make_pair(2, 3));

    ttm::smear_ttm4x smr; // smearing for chloride?
    smr.m_aDD_intra_12 = 0.626;
    smr.m_aDD_intra_13 = 0.055;

    m_electrostatics(system_nsites, system_charge, system_polfac, system_atmpolar, system_sitecrds, m_excluded12, m_excluded13, smr, Eelec, Eind, grd);

        // free up memory
        delete[] system_sitecrds;
        delete[] system_charge;
        delete[] system_polfac;
        delete[] system_atmpolar;
        delete[] h2o_atmcrds;
        delete[] x_atmcrds;

//    for (int i = 0; i < 18; ++i)
//        grd[i] = gEelec[i] + gEind[i];

#ifdef DEBUG_ENERGY 
	std::cout << "Eelec = " << Eelec << std::endl;
	std::cout << "Eind = " << Eind << std::endl;
#endif
	
	return Eelec + Eind + dispersion(crd, crd + 9, grd, grd + 9);
}

//----------------------------------------------------------------------------//

double x2b_h2o_ion_base::dispersion(const double* w1, const double*  x) const
{
    const double* O  = w1;
    const double* H1 = w1 + 3;
    const double* H2 = w1 + 6;

    const double* X  =  x;

    const double XH68 =
        x68(C6_XH, m_d6_XH, C8_XH, m_d8_XH, X, H1)
	+ x68(C6_XH, m_d6_XH, C8_XH, m_d8_XH, X, H2);

    const double XO68 =
	x68(C6_XO, m_d6_XO, C8_XO, m_d8_XO, X, O);

    return XH68 + XO68;
}

//----------------------------------------------------------------------------//

double x2b_h2o_ion_base::dispersion
    (const double* w1, const double*  x, double* g1, double* g2) const
{
    const double* O  = w1;
    const double* H1 = w1 + 3;
    const double* H2 = w1 + 6;

    const double* X  =  x;

    double* gO  = g1;
    double* gH1 = g1 + 3;
    double* gH2 = g1 + 6;

    double* gX  = g2;

    const double XH68 =
        x68(C6_XH, m_d6_XH, C8_XH, m_d8_XH, X, H1, gX, gH1)
      + x68(C6_XH, m_d6_XH, C8_XH, m_d8_XH, X, H2, gX, gH2);

    const double XO68 =
        x68(C6_XO, m_d6_XO, C8_XO, m_d8_XO, X, O, gX, gO);

    return XH68 + XO68;
}

////////////////////////////////////////////////////////////////////////////////

} // namespace h2o_ion

////////////////////////////////////////////////////////////////////////////////
