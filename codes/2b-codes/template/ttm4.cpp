#include "ttm4.h"

#include <iostream>

//TODO: merge constants.h and fit-constants.h appropriately
#include "fit-constants.h"

namespace {

const double gammaM = 0.426706882;
const double gamma1 = 1.0 - gammaM;
const double gamma2 = gammaM/2;

const double CHARGECON = constants::CHARGECON;

inline void compute_M_site_crd
    (const double O[3], const double H1[3], const double H2[3], double M[3])
{
    for (size_t i = 0; i < 3; ++i)
        M[i] = gamma1*O[i] + gamma2*(H1[i] + H2[i]);
}

} // namespace

namespace x2o {

ttm4::ttm4() {

}

ttm4::~ttm4() {
	delete[] memory;
}

ttm4::ttm4(double* crd) {
    nsites = 4;
    allocate();
    sitecrds = set_sitecrds(crd);
    atmpolar = set_pol();
    charge = set_charges(crd);
    polfac = set_polfacs(atmpolar);

    excluded12.clear();
    excluded13.clear();

    excluded12.insert(std::make_pair(0, 1)); // O - H1
    excluded12.insert(std::make_pair(0, 2)); // O - H2
    excluded12.insert(std::make_pair(0, 3)); // O - M
    excluded13.insert(std::make_pair(1, 2)); // H1 - H2
    excluded12.insert(std::make_pair(1, 3)); // H1 - M
    excluded12.insert(std::make_pair(2, 3)); // H2 - M
}

double* ttm4::set_sitecrds(double* atmcrds) {
    sitecrds = memory + nsites;
    // assumes O H H 
    compute_M_site_crd(atmcrds, atmcrds + 3, atmcrds + 6, sitecrds + 9);
    std::copy(atmcrds, atmcrds + 9, sitecrds);
    return sitecrds;
}

double* ttm4::set_charges(double* atmcrds) {
    double chgtmp[3];
    h2o::ps::dms_nasa(0.0, 0.0, 0.0, atmcrds, chgtmp, 0, false);
    const double tmp = 0.5*gammaM/(1.0 - gammaM);
    charge = memory; 

    charge[0] = 0.0;                        // O
    charge[1] = CHARGECON*(chgtmp[1] + tmp*(chgtmp[1] + chgtmp[2])); // H1
    charge[2] = CHARGECON*(chgtmp[2] + tmp*(chgtmp[1] + chgtmp[2])); // H2
    charge[3] = CHARGECON*(chgtmp[0]/(1.0 - gammaM));       // M

    //std::cout << "charge[1] = " << charge[1] << std::endl;

    return charge;
}

double* ttm4::set_pol() {
    atmpolar = memory + nsites + nsites*3;
    atmpolar[0] = 1.310; // polarO
    atmpolar[1] = 0.294; // polarH
    atmpolar[2] = 0.294; // polarH
    atmpolar[3] = 0.0;   // polarM
    return atmpolar;
}

double* ttm4::set_polfacs(double* atmpol) {
    polfac = memory + nsites + nsites*3 + nsites;
    polfac[0] = 1.310; // polarO
    polfac[1] = 0.294; // polarH
    polfac[2] = 0.294; // polarH
    polfac[3] = 1.31;   // polarM
    return polfac;
}

void ttm4::allocate() {
    memory = new double [nsites // charges
	+ nsites*3              // sitecrds	
	+ nsites                // polarizabilities
	+ nsites];              // polfacs
}

int ttm4::get_nsites() { return nsites; }
double* ttm4::get_sitecrds() { return sitecrds; }
double* ttm4::get_charges() { return charge; }
double* ttm4::get_polfacs() { return polfac; }
double* ttm4::get_pol() { return atmpolar; }

excluded_set_type::iterator ttm4::get_begin_12() { return excluded12.begin(); }
excluded_set_type::iterator ttm4::get_begin_13() { return excluded13.begin(); }
excluded_set_type::iterator ttm4::get_begin_14() { return excluded14.begin(); }
excluded_set_type::iterator ttm4::get_end_12() { return excluded12.end(); }
excluded_set_type::iterator ttm4::get_end_13() { return excluded13.end(); }
excluded_set_type::iterator ttm4::get_end_14() { return excluded14.end(); }


} // namespace x2o
