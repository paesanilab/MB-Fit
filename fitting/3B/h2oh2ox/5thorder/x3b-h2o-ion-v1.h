#ifndef X3B_H2O_ION_V2_H
#define X3B_H2O_ION_V2_H

#include <string>
#include <iostream>

#include "x3b-h2o-ion-base.h"

namespace h2o_ion {

struct x3b_h2o_ion_v1 : public x3b_h2o_ion_base {

    x3b_h2o_ion_v1();
    ~x3b_h2o_ion_v1();

    static const unsigned num_nonlinear_params = 14; 

    void set_nonlinear_parameters(const double*);
    void get_nonlinear_parameters(double*) const;

    bool nonlinear_parameters_out_of_range() const;

    static std::string name();

    // fitting interface
    static size_t nparams(); // linear parameters
    double basis(const double xyz[21], double*) const;

    // value
    void set(const double*);
    double value(const double xyz[21]) const;

    // IO
    void as_cdl(std::ostream&) const;
    void load_netcdf(const char*);

    void set_r3i_r3f(const double&, const double&);

public:
    double f_switch_cos(const double&) const;
    double f_switch_cos(const double x[21]) const;
    void doevpoly(double* x, double* p);

private:
    double m_kOH_intra;
    double m_kHH_intra;

    double m_kOO;
    double m_kOH;
    double m_kHH;

    double m_kXO; 
    double m_kXH;

    double m_dOH_intra;
    double m_dHH_intra;

    double m_dOO;
    double m_dOH;
    double m_dHH;

    double m_dXO; 
    double m_dXH;

private:
    double* m_poly;
};

} // namespace h2o_ion

#endif // X3B_H2O_ION_V2_H
