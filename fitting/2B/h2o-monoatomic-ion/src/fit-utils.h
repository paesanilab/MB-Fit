#ifndef FIT_UTILS_H
#define FIT_UTILS_H

#include "training-set-h2o-ion.h"

namespace h2o_ion {

void setup_weights(const std::vector<h2o_ion::dimer>&,
                   const double& E_range, double& E_min,
                   double* weights, size_t& effective_size,
                   bool use_mb_energy = true);

} // namespace h2o_ion

#endif // FIT_UTILS_H
