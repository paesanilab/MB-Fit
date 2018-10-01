#ifndef POLY_MODEL_H
#define POLY_MODEL_H

namespace mb_system {

struct poly_model {
    static const unsigned n_vars = 21;
    static const unsigned size = 492;

    static double eval(const double a[492],
                       const double x[21]);

    static double eval(const double a[492],
                       const double x[21],
                             double g[21]);

    static double eval_direct(const double a[492],
                              const double x[21]);

public:
    unsigned report_nvars(){ return n_vars; };
    unsigned report_size(){ return size; };
};

} // namespace mb_system

#endif // POLY_MODEL_H