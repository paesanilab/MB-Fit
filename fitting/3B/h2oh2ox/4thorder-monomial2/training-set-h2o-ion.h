#ifndef TRAINING_SET_H
#define TRAINING_SET_H

#include <cstdlib>
#include <vector>

namespace h2o_ion {

struct dimer {
    double energy_total;
    double energy_twobody;
    double energy_onebody[2];
    double xyz[12]; // O H H X 
    double mb_total_energy() const
    {
        return energy_twobody + energy_onebody[0] + energy_onebody[1];
    }
};

size_t load_dimers(const char* filename, std::vector<dimer>& ts);

struct trimer {
    double energy_total;
    double energy_threebody;
    double xyz[21]; // O H H O H H X 
};

size_t load_trimers(const char* filename, std::vector<trimer>& ts);

} // namespace h2o_ion

#endif // TRAINING_SET_H
