
#ifndef POLY_2B_A1B2Z2_D1E2_H
#define POLY_2B_A1B2Z2_D1E2_H

namespace mb_system_fit {

struct poly_model{
    static const unsigned degree = 3;
    static const unsigned n_vars = 21;

    static const unsigned size = 492;

    void eval(const double x[n_vars], double*);
};

} // namespace

#endif

