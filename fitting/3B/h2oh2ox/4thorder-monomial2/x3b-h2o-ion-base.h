#ifndef X3B_H2O_ION_BASE_H
#define X3B_H2O_ION_BASE_H

#include <iostream>
#include <vector>

#include "electrostatics.h"
#include "ttm4.h"
#include "ion.h"
#include "coulomb.h"

namespace h2o_ion {

struct x3b_h2o_ion_base {

    x3b_h2o_ion_base();

    const double& r3i() const;
    const double& r3f() const;

    double f_switch(const double&) const;
    double f_switch(const double xyz[21]) const;

    double Eind3(const double xyz[21]) const;

    virtual void as_cdl(std::ostream&) const;
    virtual void from_cdf(int ncid);

private:
    mutable ttm::electrostatics m_electrostatics;

protected:
    double m_r3i;
    double m_r3f;
};

} // namespace h2o_ion

#endif // X3B_H2O_ION_BASE_H
