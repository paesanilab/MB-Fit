#ifndef FIT_UTILS_H
#define FIT_UTILS_H

#include "training_set.h"

namespace tset {

void setup_weights(const std::vector<tset::monomer>&,
                   const double& E_range, double& E_min,
                   double* weights, size_t& effective_size,
                   bool use_mb_energy = true);

} // namespace tset

#endif // FIT_UTILS_H
