#ifndef POLY_MODEL_H
#define POLY_MODEL_H

namespace mb_system {

struct poly_model {
    static const unsigned n_vars = 31;
    static const unsigned size = 1153;

    static double eval(const double a[1153],
                       const double x[31]);

    static double eval(const double a[1153],
                       const double x[31],
                             double g[31]);

    static double eval_direct(const double a[1153],
                              const double x[31]);

public:
    unsigned report_nvars(){ return n_vars; };
    unsigned report_size(){ return size; };
};

} // namespace mb_system

#endif // POLY_MODEL_H