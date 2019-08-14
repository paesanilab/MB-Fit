#ifndef POLY_MODEL_H
#define POLY_MODEL_H

namespace mb_system {

struct poly_model {
    static const unsigned n_vars = 3;
    static const unsigned size = 6;

    static double eval(const double a[6],
                       const double x[3]);

    static double eval(const double a[6],
                       const double x[3],
                             double g[3]);

    static double eval_direct(const double a[6],
                              const double x[3]);

public:
    unsigned report_nvars(){ return n_vars; };
    unsigned report_size(){ return size; };
};

} // namespace mb_system

#endif // POLY_MODEL_H