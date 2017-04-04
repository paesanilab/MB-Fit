#ifndef POLY_MODEL_H
#define POLY_MODEL_H

namespace mb_system {

struct poly_model {
    static const unsigned n_vars = 45;
    static const unsigned size = 53;

    static double eval(const double a[53],
                       const double x[45]);

    static double eval(const double a[53],
                       const double x[45],
                             double g[45]);

    static double eval_direct(const double a[53],
                              const double x[45]);

public:
    unsigned report_nvars(){ return n_vars; };
    unsigned report_size(){ return size; };
};

} // namespace mb_system

#endif // POLY_MODEL_H