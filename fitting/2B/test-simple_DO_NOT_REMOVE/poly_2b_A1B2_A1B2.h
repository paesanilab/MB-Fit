#ifndef POLY_2B_A1B2_A1B2_H
#define POLY_2B_A1B2_A1B2_H

namespace mb_system_fit {

struct poly_model{ 
    static const unsigned degree = 3;
    static const unsigned n_vars = 15;

    static const unsigned size = 134;

    void eval(const double x[n_vars], double*);
};

} // namespace

#endif
