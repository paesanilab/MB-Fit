#ifndef X1B_V1_H
#define X1B_V1_H

#include <string>
#include <iostream>

namespace h3o {

struct x1b_v1 {

    x1b_v1();
    ~x1b_v1();

    static const unsigned num_nonlinear_params = 4;

    void set_nonlinear_parameters(const double*);
    void get_nonlinear_parameters(double*) const;

    bool nonlinear_parameters_out_of_range() const;

    static std::string name();

    // fitting interface
    static size_t nparams(); // linear parameters
    double basis(const double xyz[12], double*) const;

    // value
    void set(const double*);
    double value(const double xyz[12]) const;

    // IO
    void as_poly_dat(std::ostream&) const;
    void load_poly_dat(const char*);

//    void set_r3i_r3f(const double&, const double&);

public:
//    double f_switch_cos(const double&) const;
//    double f_switch_cos(const double x[27]) const;

private:
    double m_kOH_intra;
    double m_kHH_intra;

    double m_dOH_intra;
    double m_dHH_intra;

private:
    double* m_poly;
};

} // namespace h3o

#endif // X1B_V1_H
