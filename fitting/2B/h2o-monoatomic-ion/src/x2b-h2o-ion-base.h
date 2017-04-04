#ifndef H2O_ION_BASE_H
#define H2O_ION_BASE_H

#include <map>
#include <string>
#include <iostream>

#include "electrostatics.h"
#include "ttm4.h"
#include "ion.h"
#include "coulomb.h"
#include <vector>

namespace h2o_ion { 

//
//  ttm4-es + TT6*C6/R^6 + TT8*C8/R^8
//

struct x2b_h2o_ion_base {

    x2b_h2o_ion_base();
    virtual ~x2b_h2o_ion_base() {} //virtual members can be modified in derived classes

    static const double C6_XH;
    static const double C6_XO;

    static const double C8_XH;
    static const double C8_XO;
    
/*    static const double C6_HH; 
    static const double C6_OH;
    static const double C6_OO;

    static const double C8_HH;
    static const double C8_OH;
    static const double C8_OO;
  */  

    void setup(const char*);

    // I/O
    virtual void as_cdl(std::ostream&) const;
    virtual void from_cdf(int ncid);

    // main functionality
    double value(const double*) const; // dispersion() + Eind + Eelec
    double dispersion(const double*) const;

    // not necessary for fitting
    // if you decide to use this, FIXME: has to set Edisp to 0
    void components(const double*,
        double& Eelec, double& Eind, double& Edisp) const;

    // including the gradients
    double value(const double*, double*) const;
    double dispersion(const double*, const double*) const;
    double dispersion(const double*, const double*, double*, double*) const;

protected:
    virtual void do_setup(const std::map<std::string, double>&);

protected:
    double m_d6_XH;
    double m_d6_XO;

    double m_d8_XH;
    double m_d8_XO;

    
/*protected:
    double m_d6_HH;
    double m_d6_OH;
    double m_d6_OO;

    double m_d8_HH;
    double m_d8_OH;
    double m_d8_OO;
*/
protected:
    mutable ttm::electrostatics m_electrostatics;
};
} // namespace h2o_ion

#endif // H2O_ION_BASE_H
