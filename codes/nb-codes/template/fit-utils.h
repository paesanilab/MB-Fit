#ifndef FIT_UTILS_H
#define FIT_UTILS_H

#include "training_set.h"

namespace tset {

void setup_weights(const std::vector<tset::nb_system>&,
                   const double& E_range, double& E_min, size_t &index_min,
                   std::vector<double> &weights, size_t& effective_size);

} // namespace tset

#endif // FIT_UTILS_H
