#ifndef FIT_UTILS_H
#define FIT_UTILS_H

#include "training-set.h"

namespace h3o {

void setup_weights(const std::vector<h3o::monomer>&,
                   const double& E_range, double& E_min,
		   double* weights, size_t& effective_size);

} // namespace tools

#endif // FIT_UTILS_H
