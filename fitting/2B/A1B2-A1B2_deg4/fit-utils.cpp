#include <cmath>
#include <cassert>

#include "fit-utils.h"

namespace tset {

void setup_weights(const std::vector<tset::dimer>& ts,
                   const double& E_range, double& E_min,
                   double* weights, size_t& effective_size,
                   bool use_mb_energy)
{
    E_min = 1.0e+10;

    for (size_t n = 0; n < ts.size(); ++n) {
        const double energy =
            (use_mb_energy ? ts[n].mb_total_energy() : ts[n].energy_total);
        if (energy < E_min)
            E_min = energy;
    }

    for (size_t n = 0; n < ts.size(); ++n) {
        const double energy =
            (use_mb_energy ? ts[n].mb_total_energy() : ts[n].energy_total);
        weights[n] = std::pow(E_range/(energy - E_min + E_range), 2);
    }

    double weights_sum(0);
    for (size_t n = 0; n < ts.size(); ++n)
        weights_sum += weights[n];

    weights_sum /= ts.size();
    double weights_sq_sum(0);
    for (size_t n = 0; n < ts.size(); ++n) {
        const double ww = weights[n]/weights_sum;
        weights[n] = ww;
        weights_sq_sum += ww*ww;
    }

    weights_sq_sum /= ts.size();
    weights_sq_sum /= ts.size();

    effective_size = 1.0/weights_sq_sum;
}

} // namespace tset
