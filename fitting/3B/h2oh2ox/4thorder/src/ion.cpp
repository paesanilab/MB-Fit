#include "ion.h"
#include "ps.h"
#include <algorithm>
#include <cstddef>
#include "fit-constants.h"

#define FLUORIDE    
//#define CHLORIDE    
//#define BROMIDE     
//#define IODIDE      

const double CHARGECON = constants::CHARGECON;

namespace ion {

    ion::ion() {

    }

    ion::~ion() {
		delete[] memory;
    }

    ion::ion( double* crd) {

	nsites = 1;
	allocate();
	sitecrds = set_sitecrds(crd);
	atmpolar = set_pol();
	charge = set_charges(crd);
	polfac = set_polfacs(atmpolar);

	//excluded pairs stuff - not applicable
	
    }

    double* ion::set_charges(double* atmcrds) {
	charge = memory;
	charge[0] = -1.0*CHARGECON;
	return charge;
    }

    double* ion::set_pol() {
	atmpolar = memory + nsites + nsites*3;
        
#ifdef FLUORIDE     
	atmpolar[0] = 2.4669; //A^3 from finite field method 
#endif

#ifdef CHLORIDE     
	atmpolar[0] = 5.3602; //A^3 from finite field method 
#endif

#ifdef BROMIDE      
	atmpolar[0] = 7.1668; //A^3 from finite field method 
#endif

#ifdef IODIDE       
	atmpolar[0] =10.1184; //A^3 from finite field method 
#endif

	return atmpolar;
    }

    double* ion::set_polfacs(double* atmpol) {
	polfac = memory + nsites + nsites*3 + nsites;
	std::copy(atmpolar, atmpolar + nsites, polfac); //polfac set equal to polar
	return polfac;
    }

    double* ion::set_sitecrds(double* atmcrds) {
	sitecrds = memory + nsites;
	return atmcrds; //no additional sites, just ion 
    }

    void ion::allocate() {

	memory = new double [nsites  //charge
	    + nsites*3  //site coordinates
	    + nsites  //polarizabilities 
	    + nsites];  //polfac
    }

    int ion::get_nsites() { return nsites; }
    double* ion::get_charges() { return charge; }
    double* ion::get_sitecrds() { return sitecrds; }
    double* ion::get_pol() { return atmpolar; }
    double* ion::get_polfacs() { return polfac; }

    excluded_set_type::iterator ion::get_begin_12() {  }
    excluded_set_type::iterator ion::get_begin_13() {  }
    excluded_set_type::iterator ion::get_end_12() {  }
    excluded_set_type::iterator ion::get_end_13() {  }

} // namespace ion

