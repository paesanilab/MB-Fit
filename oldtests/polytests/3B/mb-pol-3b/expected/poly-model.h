#ifndef POLY_MODEL_H
#define POLY_MODEL_H

namespace mb_system {

struct poly_model {
    static const unsigned n_vars = 36;
    static const unsigned size = 1163;

    static double eval(const double a[1163],
                       const double x[36]);

    static double eval(const double a[1163],
                       const double x[36],
                             double g[36]);

    static double eval_direct(const double a[1163],
                              const double x[36]);

public:
    unsigned report_nvars(){ return n_vars; };
    unsigned report_size(){ return size; };
};

} // namespace mb_system

#endif // POLY_MODEL_H