#ifndef POLY_2B_H2O_ION_H
#define POLY_2B_H2O_ION_H

namespace h2o_ion {

//
// reduced (no bonds breaking) permutational symmetry
// including only 2B terms
//
//    x[0] = @VAR@-|x-intra|(H1, H2);
//    x[1] = @VAR@-|x-intra|(O, H1);
//    x[2] = @VAR@-|x-intra|(O, H2);
//    x[3] = @VAR@-|x-coul|(X, H1);
//    x[4] = @VAR@-|x-coul|(X, H2);
//    x[5] = @VAR@-|x-coul|(X, O);
//    x[6] = @VAR@-|x-main|(X, L1);
//    x[7] = @VAR@-|x-main|(X, L2);
// specialized for D = 1/2/3 (what does this mean?) 
//

template <unsigned D>
struct poly_2b_h2o_ion {
    static const unsigned degree = D;
    static const unsigned n_vars = 8; //8 ksis for water-ion system
    static const unsigned size;

    static void eval(const double x[n_vars], double*);
};

} // namespace

#endif // POLY_2B_H2O_ION_H
