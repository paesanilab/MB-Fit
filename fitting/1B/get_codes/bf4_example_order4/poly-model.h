#ifndef POLY_MODEL_H
#define POLY_MODEL_H

namespace mb_system {

struct poly_model {
    static const unsigned n_vars = 10;
    static const unsigned size = 82;

    static double eval(const double a[82],
                       const double x[10]);

    static double eval(const double a[82],
                       const double x[10],
                             double g[10]);

    static double eval_direct(const double a[82],
                              const double x[10]);

public:
    unsigned report_nvars(){ return n_vars; };
    unsigned report_size(){ return size; };
};

} // namespace mb_system

#endif // POLY_MODEL_H