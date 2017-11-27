#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif // HAVE_CONFIG_H

#include <cassert>
#include <cstdlib>

#include <netcdf.h>

#include <iomanip>
#include <iostream>
#include <algorithm>

#include "stuff.h"    //what's in stuff?
#include "x3b-h2o-ion-base.h"

/////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

void error(int ec)
{
    std::cerr << " ** Fatal Error in x3b_h2o_ion_base::from_cdf() ** : "
              << nc_strerror(ec) << std::endl;
    std::exit(EXIT_FAILURE);
}

//-----------------------------------------------------------------------------//

} // namespace

/////////////////////////////////////////////////////////////////////////////////

namespace h2o_ion {

//-----------------------------------------------------------------------------//

x3b_h2o_ion_base::x3b_h2o_ion_base()
{
    //set appropriate values for the switching distances

    m_r3i = 0.0;          
    m_r3f = 6.5;
}

//-----------------------------------------------------------------------------//

const double& x3b_h2o_ion_base::r3i() const
{
    return m_r3i;
}

//-----------------------------------------------------------------------------//

const double& x3b_h2o_ion_base::r3f() const
{
    return m_r3f;
}

//-----------------------------------------------------------------------------//

double x3b_h2o_ion_base::f_switch(const double& r) const
{
    if (r > m_r3f) {
        return 0.0;
    } else if (r > m_r3i) {
        const double t1 = 1.0/(m_r3f - m_r3i);
        const double x = (r - m_r3i)*t1;
        return 1 + x*x*(2*x - 3.0);
    } else {
        return 1.0;
    }
}

//-----------------------------------------------------------------------------//

double x3b_h2o_ion_base::f_switch(const double xyz[21]) const
{
    const double s12 = f_switch(distance(xyz + 0, xyz + 9));
    const double s13 = f_switch(distance(xyz + 0, xyz + 18));
    const double s23 = f_switch(distance(xyz + 9, xyz + 18));

    const double s = s12*s13 + s12*s23 + s13*s23;

    return s;
}

//-----------------------------------------------------------------------------//

double x3b_h2o_ion_base::Eind3(const double xyz[21]) const
{
    //The electrostatics stuff.
    double* h2o1_atmcrds = new double [9];
    std::copy(xyz, (xyz + 9), h2o1_atmcrds);

    double* h2o2_atmcrds = new double [9];
    std::copy((xyz + 9), (xyz + 18), h2o2_atmcrds);

    double* x_atmcrds = new double [3];
    std::copy((xyz + 18), (xyz + 21), x_atmcrds);


    x2o::ttm4 h2o1(h2o1_atmcrds);
    x2o::ttm4 h2o2(h2o2_atmcrds);
    ion::ion x1(x_atmcrds);
    std::vector<molecule*> system;
    system.push_back(&h2o1);
    system.push_back(&h2o2);
    system.push_back(&x1);

    int system_nsites = 0;
    double* system_sitecrds;
    double* system_charge;
    double* system_polfac;
    double* system_atmpolar;

    for (int i = 0; i < system.size(); i++) {
      system_nsites += system[i]->get_nsites();
    } //calculates total number of sites in the system

    //allocate memory to the pointers
    system_sitecrds = new double [system_nsites*3];
    system_charge = new double [system_nsites];
    system_polfac = new double [system_nsites];
    system_atmpolar = new double [system_nsites]; 

    int pos = 0;
    for(int i = 0; i < system.size(); i++) {
        
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
    m_excluded12.insert(std::make_pair(4, 5));
    m_excluded12.insert(std::make_pair(4, 6));
    m_excluded12.insert(std::make_pair(4, 7));
    m_excluded13.insert(std::make_pair(1, 2));
    m_excluded13.insert(std::make_pair(5, 6));
    m_excluded12.insert(std::make_pair(1, 3));
    m_excluded12.insert(std::make_pair(2, 3));
    m_excluded12.insert(std::make_pair(5, 7));
    m_excluded12.insert(std::make_pair(6, 7));
    //write a snippet to make these excluded pairs for a general system at some point

    ttm::smear_ttm4x smr; 
    smr.m_aDD_intra_12 = 0.626;
    smr.m_aDD_intra_13 = 0.055;

    double Eind12, Eind13, Eind23, Eelec, Eind;

    m_electrostatics(9, system_charge, system_polfac, system_atmpolar, system_sitecrds, m_excluded12, m_excluded13, smr, Eelec, Eind, 0);

    m_electrostatics(8, system_charge, system_polfac, system_atmpolar, system_sitecrds, m_excluded12, m_excluded13, smr, Eelec, Eind12, 0);

    m_electrostatics(5, (system_charge + 4), (system_polfac + 4), (system_atmpolar + 4), (system_sitecrds + 12), m_excluded12, m_excluded13, smr, Eelec, Eind23, 0); 

    double* sitecrds13 = new double [15];
    std::copy(system_sitecrds, (system_sitecrds + 12), sitecrds13);
    std::copy((system_sitecrds + 24), (system_sitecrds + 27), (sitecrds13 + 12));
    
    double* charge13 = new double [5];
    std::copy(system_charge, (system_charge + 4), charge13);
    std::copy((system_charge + 8), (system_charge + 9), (charge13 + 4));
    
    double* polfac13 = new double [5];
    std::copy(system_polfac, (system_polfac + 4), polfac13);
    std::copy((system_polfac + 8), (system_polfac + 9), (polfac13 + 4));
    
    double* atmpolar13 = new double [5];
    std::copy(system_atmpolar, (system_atmpolar + 4), atmpolar13);
    std::copy((system_atmpolar + 8), (system_atmpolar + 9), (atmpolar13 + 4));
    
    m_electrostatics(5, charge13, polfac13, atmpolar13, sitecrds13, m_excluded12, m_excluded13, smr, Eelec, Eind13, 0);
    /*double Eind12, Eind13, Eind23, Eelec, Eind;

    m_ttm4(3, xyz, Eelec, 0, Eind, 0);

    m_ttm4(2, xyz, Eelec, 0, Eind12, 0);
    m_ttm4(2, xyz + 9, Eelec, 0, Eind23, 0);

    double xyz13[18];

    std::copy(xyz, xyz + 9, xyz13);
    std::copy(xyz + 18, xyz + 27, xyz13 + 9);

    m_ttm4(2, xyz13, Eelec, 0, Eind13, 0);
    */

    //free up memory
    delete[] system_sitecrds;
    delete[] system_charge;
    delete[] system_polfac;
    delete[] system_atmpolar;
    delete[] h2o1_atmcrds;
    delete[] h2o2_atmcrds;
    delete[] x_atmcrds;
    delete[] sitecrds13;
    delete[] charge13;
    delete[] polfac13;
    delete[] atmpolar13;


    return Eind - Eind12 - Eind23 - Eind13; 
 
}

//-----------------------------------------------------------------------------//

void x3b_h2o_ion_base::as_cdl(std::ostream& os) const
{
    using namespace std;

    ios::fmtflags saved_flags = os.flags();

    os << setprecision(15) << scientific
       << "  :r3i = " << setw(22) << m_r3i << "; // A\n"
          "  :r3f = " << setw(22) << m_r3f << "; // A\n";

    os.flags(saved_flags);
}

//-----------------------------------------------------------------------------//

void x3b_h2o_ion_base::from_cdf(int ncid)
{
    int rc;

#   define RETRIEVE(name) \
    if ((rc = nc_get_att_double(ncid, NC_GLOBAL, #name , &m_##name))) \
        error(rc);

    RETRIEVE(r3i)
    RETRIEVE(r3f)

#   undef RETRIEVE
}

//-----------------------------------------------------------------------------//

} // namespace h2o_ion

/////////////////////////////////////////////////////////////////////////////////
