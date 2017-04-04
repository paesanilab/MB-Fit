#ifndef POLY_MODEL_H
#define POLY_MODEL_H

namespace mb_system {

struct poly_model {
    static const unsigned n_vars = 6;
    static const unsigned size = 50;

    static double eval(const double a[50],
                       const double x[6]);

    static double eval(const double a[50],
                       const double x[6],
                             double g[6]);

    static double eval_direct(const double a[50],
                              const double x[6]);

public:
    unsigned report_nvars(){ return n_vars; };
    unsigned report_size(){ return size; };
};

} // namespace mb_system

#endif // POLY_MODEL_H