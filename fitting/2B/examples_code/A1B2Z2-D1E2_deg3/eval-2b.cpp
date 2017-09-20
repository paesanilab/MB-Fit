
#include <cmath>
#include <cassert>
#include <cstdlib>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "mon1.h"
#include "mon2.h"
#include "training_set.h"
#include "x2b_A1B2Z2_D1E2_v1x.h"
#include "electrostatics.h"
#include "coulomb.h"
#include "dispersion.h"
#include "io-xyz.h"

#define GRADIENTS

static std::vector<double> elec_e;
static std::vector<double> disp_e;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: eval fit-2b.nc dimer.xyz"
                  << std::endl;
        return 0;
    }
    std::cout << std::scientific << std::setprecision(9);
    x2b_A1B2Z2_D1E2::x2b_A1B2Z2_D1E2_v1x pot;
    std::vector<std::string> elements;
    std::vector<double> crd;

    try {
        ++argv;
        --argc;
        pot.load_netcdf(*argv);

        ++argv;
        --argc;
        std::ifstream ifs(*argv);

        if(!ifs)
            throw std::runtime_error("could not open the XYZ file");

        std::string comment;
        kit::io::load_xyz(ifs, comment, elements, crd);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }
    
    double xyz[18];
    std::copy(crd.begin(), crd.end(), xyz);
    
    x::mon1 m1(xyz);
    x::mon2 m2(xyz + 3*m1.get_realsites());
    
    int system_nsites = m1.get_nsites() + m2.get_nsites();
    int * system_is_w;
    double* system_sitecrds;
    double* system_charge;
    double* system_polfac;
    double* system_pol;

    system_sitecrds = new double [system_nsites*3];
    system_charge = new double [system_nsites];
    system_polfac = new double [system_nsites];
    system_pol = new double [system_nsites]; //allocates memory to the pointers
    system_is_w = new int[system_nsites];
      
    std::fill(system_is_w, system_is_w + m1.get_nsites(), m1.is_w);
    std::fill(system_is_w + m1.get_nsites(), system_is_w + system_nsites, m2.is_w);

    std::copy(m1.get_sitecrds(), m1.get_sitecrds() + 3 * m1.get_nsites(),
              system_sitecrds);
    std::copy(m2.get_sitecrds(), m2.get_sitecrds() + 3 * m2.get_nsites(),
              system_sitecrds + 3 * m1.get_nsites());

    std::copy(m1.get_charges(), m1.get_charges() + m1.get_nsites(),
              system_charge);
    std::copy(m2.get_charges(), m2.get_charges() + m2.get_nsites(),
              system_charge + m1.get_nsites());

    std::copy(m1.get_pol(), m1.get_pol() + m1.get_nsites(),
              system_pol);
    std::copy(m2.get_pol(), m2.get_pol() + m2.get_nsites(),
              system_pol + m1.get_nsites());

    std::copy(m1.get_polfacs(), m1.get_polfacs() + m1.get_nsites(),
              system_polfac);
    std::copy(m2.get_polfacs(), m2.get_polfacs() + m2.get_nsites(),
              system_polfac + m1.get_nsites());

    excluded_set_type exclude12;
    excluded_set_type exclude13;
    excluded_set_type exclude14;

    for (auto i = m1.get_begin_12(); i != m1.get_end_12(); i++) {
        exclude12.insert(*i);
    }
    for (auto i = m2.get_begin_12(); i != m2.get_end_12(); i++) {
        std::pair<size_t,size_t> p =
                    std::make_pair(i->first + m1.get_nsites() ,
                                   i->second + m1.get_nsites());
        exclude12.insert(p);
    }

    for (auto i = m1.get_begin_13(); i != m1.get_end_13(); i++) {
        exclude13.insert(*i);
    }
    for (auto i = m2.get_begin_13(); i != m2.get_end_13(); i++) {
      std::pair<size_t,size_t> p =
                    std::make_pair(i->first + m1.get_nsites() ,
                                   i->second + m1.get_nsites());
        exclude13.insert(p);
    }

    for (auto i = m1.get_begin_14(); i != m1.get_end_14(); i++) {
        exclude14.insert(*i);
    }
    for (auto i = m2.get_begin_14(); i != m2.get_end_14(); i++) {
      std::pair<size_t,size_t> p =
                    std::make_pair(i->first + m1.get_nsites() ,
                                   i->second + m1.get_nsites());
        exclude14.insert(p);
    }

    ttm::electrostatics m_electrostatics;
    ttm::smear_ttm4x smr; 
    smr.m_aDD_intra_12 = 0.3;
    smr.m_aDD_intra_13 = 0.3;
    smr.m_aDD_intra_14 = 0.055;

    double ener = m_electrostatics(system_nsites, system_charge, system_polfac, system_pol,
                                   system_sitecrds, exclude12, exclude13, exclude14, 
                                   system_is_w, smr, 0);

    elec_e.push_back(ener);
      
    // Now need to take out dispersion
    x2b_disp disp(m1.get_sitecrds(), m2.get_sitecrds(), m1.get_realsites(), m2.get_realsites());
    ener = disp.get_dispersion();
    disp_e.push_back(ener);
    
    double Epoly = pot(xyz);
    std::cout << "E_nograd = " << Epoly + elec_e[0] + disp_e[0] << std::endl;
    std::cout << "E_poly = " << Epoly << std::endl;
    std::cout << "E_elec = " << elec_e[0] << std::endl;
    std::cout << "E_disp = " << disp_e[0] << std::endl;
    
#ifdef GRADIENTS
    const double eps = 1.0e-5;
    double grd[18];
    std::fill(grd, grd + 18, 0.0);
    Epoly = pot(xyz, grd);
    std::cout << "E_grad = " << Epoly << std::endl;
    for(size_t n = 0; n < 18; ++n){
      const double x_orig = xyz[n];

      xyz[n] = x_orig + eps;
      const double Ep = pot(xyz) ;
      std::cout << "E_p = " << Ep << std::endl;

      xyz[n] = x_orig + 2*eps;
      const double E2p = pot(xyz) ;
      std::cout << "E_2p = " << E2p << std::endl;

      xyz[n] = x_orig - 2*eps;
      const double E2m = pot(xyz) ;
      std::cout << "E_2m = " << E2m << std::endl;

      xyz[n] = x_orig - eps;
      const double Em = pot(xyz) ;
      std::cout << "E_m = " << Em << std::endl;

      const double gfd = (8*(Ep - Em) - (E2p - E2m))/(12*eps);
      xyz[n] = x_orig;

      std::cout << elements[n/3] << "   "  << "Analit: " << grd[n] << " Numerical: " << gfd
                     << " Diff: " << std::fabs(grd[n] - gfd) << '\n';
    }
#endif

    // Free memory
    delete[] system_sitecrds ;
    delete[] system_charge  ;
    delete[] system_polfac  ;
    delete[] system_pol ;
    delete[] system_is_w;
    
    return 0;
}
