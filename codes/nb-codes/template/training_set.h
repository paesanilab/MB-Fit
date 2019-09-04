#ifndef TRAINING_SET_H
#define TRAINING_SET_H

#include <cstdlib>
#include <vector>

namespace tset {

struct nb_system {
    double binding_energy;
    double nb_energy;
    
    // Note: XYZ is for the real sites
    std::vector<double> xyz;
};

size_t load_nb_system(const char* filename, std::vector<nb_system>& ts);

} // namespace tset

#endif // TRAINING_SET_H
