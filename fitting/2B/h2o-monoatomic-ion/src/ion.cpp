#include "ion.h"
#include "ps.h"
#include <algorithm>
#include <cstddef>
#include "fit-constants.h"

//#define FLUORIDE
//#define CHLORIDE
//#define BROMIDE
//#define IODIDE

//#define VERBOSE
#ifdef VERBOSE 
#include <iostream>
#endif 

const double CHARGECON = constants::CHARGECON;

namespace x  {

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
	atmpolar[0] = 2.4669; //A^3 CCSD(T)/t-av5z  from finite field method 
#ifdef POL75 
	atmpolar[0] = 1.8502; //A^3 
#endif
#ifdef POL50 
	atmpolar[0] = 1.2334; //A^3 
#endif
#ifdef POL25 
	atmpolar[0] = 0.6167; //A^3 
#endif
#ifdef POL0  
	atmpolar[0] = 0.0000; //A^3
#endif
#ifdef POL_EFF
	atmpolar[0] = 1.0651 ; //A^3 CCSD(T)/t-av5z  scaling factor 0.4318
#endif
#endif 

#ifdef CHLORIDE 
	atmpolar[0] = 5.3602; //A^3 CCSD(T)/t-av5z from finite field method 
#ifdef POL75 
	atmpolar[0] = 4.0201; //A^3 
#endif
#ifdef POL50 
	atmpolar[0] = 2.6801; //A^3 
#endif
#ifdef POL25 
	atmpolar[0] = 1.3400; //A^3 
#endif
#ifdef POL0  
	atmpolar[0] = 0.0000; //A^3
#endif
#ifdef POL_EFF
	atmpolar[0] = 3.1495 ; //A^3 CCSD(T)/t-av5z  scaling factor 0.5976 
#endif
#endif

#ifdef BROMIDE  
	atmpolar[0] = 7.1668; //A^3 CCSD(T)/t-av5z from finite field method 
#ifdef POL75 
	atmpolar[0] = 5.3751; //A^3 
#endif
#ifdef POL50 
	atmpolar[0] = 3.5834; //A^3 
#endif
#ifdef POL25 
	atmpolar[0] = 1.7917; //A^3 
#endif
#ifdef POL0  
	atmpolar[0] = 0.0000; //A^3 
#endif
#ifdef POL_EFF
	atmpolar[0] = 4.4640 ; //A^3 CCSD(T)/t-av5z  scaling factor 0.6229 
#endif
#endif

#ifdef IODIDE  
	atmpolar[0] = 10.1184; //A^3 CCSD(T)/t-av5z  from finite field method 
#ifdef POL75 
	atmpolar[0] = 7.5888; //A^3 CCSD(T)/t-av5z  from finite field method 
#endif
#ifdef POL50 
	atmpolar[0] = 5.0592; //A^3 CCSD(T)/t-av5z  from finite field method 
#endif
#ifdef POL25 
	atmpolar[0] = 2.5296; //A^3 CCSD(T)/t-av5z  from finite field method 
#endif
#ifdef POL0  
	atmpolar[0] = 0.0000; //A^3 CCSD(T)/t-av5z  from finite field method 
#endif
#ifdef POL_EFF
	atmpolar[0] = 6.9025 ; //A^3 CCSD(T)/t-av5z  scaling factor 0.6822 
#endif
#endif

#ifdef VERBOSE 
        for (size_t i=0; i<nsites; i++) {
             std::cout << "Ion polarizability : " << atmpolar[0] ; 
        } 
#endif

	return atmpolar;
    }

    double* ion::set_polfacs(double* atmpol) {
	polfac = memory + nsites + nsites*3 + nsites;
#ifdef FLUORIDE 
	polfac[0] = 2.4669;  //100% polarizability
#endif 
#ifdef CHLORIDE 
	polfac[0] = 5.3602;  //100% polarizability
#endif 
#ifdef BROMIDE  
	polfac[0] = 7.1668;  //100% polarizability
#endif 
#ifdef IODIDE   
	polfac[0] = 10.1184; //100% polarizability
#endif 
	//std::copy(atmpolar, atmpolar + nsites, polfac); //is this OK?
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

} // namespace x 





    


