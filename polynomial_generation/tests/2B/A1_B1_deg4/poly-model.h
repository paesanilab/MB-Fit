#ifndef POLY_MODEL_H
#define POLY_MODEL_H

namespace mb_system {

struct poly_model {
    static const unsigned n_vars = 1;
    static const unsigned size = 4;

    static double eval(const double a[4],
                       const double x[1]);

    static double eval(const double a[4],
                       const double x[1],
                             double g[1]);

    static double eval_direct(const double a[4],
                              const double x[1]);

public:
    unsigned report_nvars(){ return n_vars; };
    unsigned report_size(){ return size; };
};

} // namespace mb_system

#endif // POLY_MODEL_H