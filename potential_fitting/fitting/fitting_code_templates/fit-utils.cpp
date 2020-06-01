#include <cmath>
#include <cassert>

#include "fit-utils.h"

namespace tset {

void setup_weights(const std::vector<tset::nb_system>& ts,
                   const double& E_range, double& E_min, size_t &index_min,
                   std::vector<double> &weights, size_t& effective_size)
{
    E_min = 1.0e+10;
    index_min = 0;

    for (size_t n = 0; n < ts.size(); ++n) {
        const double energy = ts[n].binding_energy;
        if (energy < E_min) {
            E_min = energy;
            index_min = n;
        }
    }

    for (size_t n = 0; n < ts.size(); ++n) {
        const double energy = ts[n].binding_energy;
        weights[n] = std::pow(E_range/(energy - E_min + E_range), 2);
    }

    double weights_sum(0);
    double weights_sq_sum(0);
    for (size_t n = 0; n < ts.size(); ++n) {
        weights_sum += weights[n];
        weights_sq_sum += weights[n]*weights[n];
    }

    weights_sum /= ts.size();
    weights_sq_sum = weights_sq_sum / ts.size() / ts.size();
//    for (size_t n = 0; n < ts.size(); ++n) {
//        const double ww = weights[n]/weights_sum;
//        weights[n] = ww;
//        weights_sq_sum += ww*ww;
//    }
//
//    weights_sq_sum /= ts.size();
//    weights_sq_sum /= ts.size();

    effective_size = 1.0/weights_sq_sum;
}

} // namespace tset
