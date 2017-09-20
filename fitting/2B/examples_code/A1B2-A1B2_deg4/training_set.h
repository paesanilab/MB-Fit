
#ifndef TRAINING_SET_H
#define TRAINING_SET_H

#include <cstdlib>
#include <vector>

namespace tset {

struct dimer {
    double energy_total;
    double energy_twobody;
    double energy_onebody[2];
    
    // Note: XYZ is for the real sites
    double xyz[18];

    double mb_total_energy() const
    {
        return energy_twobody + energy_onebody[0] + energy_onebody[1];
    }
};

size_t load_dimers(const char* filename, std::vector<dimer>& ts);

} // namespace tset

#endif // TRAINING_SET_H
