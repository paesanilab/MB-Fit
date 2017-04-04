#ifndef TRAINING_SET_H
#define TRAINING_SET_H

#include <cstddef>
#include <vector>

namespace h3o {

struct monomer{
    double energy_onebody;
    double xyz[12];
};

size_t load_monomers(const char* filename, std::vector<monomer>& ts);

} // namespace x2o

#endif // TRAINING_SET_H
