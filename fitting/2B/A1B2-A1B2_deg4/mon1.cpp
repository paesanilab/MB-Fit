
#include "mon1.h"
#include "ps.h"
#include <algorithm>
#include <cstddef>
#include "fit-constants.h"

//#define VERBOSE
#ifdef VERBOSE 
#include <iostream>
#endif

const double CHARGECON = constants::CHARGECON;

// NOTE: ASSUMING ABB

namespace x  {

  mon1::mon1() { }

  mon1::~mon1() {
    delete[] memory;
  }

  mon1::mon1( double* crd) {
    
    nsites = 3; 
    realsites = 3;
    is_w = 0;
    allocate();
    
    sitecrds = set_sitecrds(crd);
    atmpolar = set_pol();
    charge = set_charges(crd);
    polfac = set_polfacs(atmpolar);

    // Excluded pairs
    excluded12.insert(std::make_pair(0,1));
    excluded12.insert(std::make_pair(0,2));
    excluded13.insert(std::make_pair(1,2));


  }

  double* mon1::set_charges(double* atmcrds) {
    charge = memory;
    charge[0] = 0.454448*CHARGECON;
    charge[1] = -0.227224*CHARGECON;
    charge[2] = -0.227224*CHARGECON;


    return charge;
  }
  
  double* mon1::set_pol() {
    atmpolar = memory + nsites + nsites*3;
    atmpolar[0] = 1.43039;
    atmpolar[1] = 0.771519;
    atmpolar[2] = 0.771519;

    return atmpolar;
  }

  double* mon1::set_polfacs(double* atmpol) {
    polfac = memory + nsites + nsites*3 + nsites;
    
    polfac[0] = 1.43039;
    polfac[1] = 0.771519;
    polfac[2] = 0.771519;

    return polfac;
  }

  double* mon1::set_sitecrds(double* atmcrds) {
    sitecrds = memory + nsites;
    
// ##DEFINE HERE## ficticious site coordinates
    return atmcrds;  
  }

  void mon1::allocate() {
    memory = new double [nsites  //charge
      + nsites*3  //site coordinates
      + nsites  //polarizabilities 
      + nsites];  //polfac
    }

    int mon1::get_nsites() { return nsites; }
    int mon1::get_realsites() { return realsites; }
    double* mon1::get_charges() { return charge; }
    double* mon1::get_sitecrds() { return sitecrds; }
    double* mon1::get_pol() { return atmpolar; }
    double* mon1::get_polfacs() { return polfac; }

    excluded_set_type::iterator mon1::get_begin_12() { return excluded12.begin(); }
    excluded_set_type::iterator mon1::get_begin_13() { return excluded13.begin(); }
    excluded_set_type::iterator mon1::get_begin_14() { return excluded14.begin(); }
    excluded_set_type::iterator mon1::get_end_12() { return excluded12.end(); }
    excluded_set_type::iterator mon1::get_end_13() { return excluded13.end(); }
    excluded_set_type::iterator mon1::get_end_14() { return excluded14.end(); }

} // namespace x
