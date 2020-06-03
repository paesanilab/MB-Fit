#include <cassert>

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "io-xyz.h"
#include "training_set.h"

namespace tset {

size_t load_nb_system(const char* filename, std::vector<nb_system>& ts, bool is_training_set)
{
    assert(filename);

    std::ifstream ifs(filename);
    if (!ifs) {
        std::ostringstream oss;
        oss << "could not open '" << filename << "' for reading";
        throw std::runtime_error(oss.str());
    }

    size_t n_nb_sys = 0;

    std::string comment;
    std::vector<std::string> elements;
    std::vector<double> xyz;

    std::istringstream iss;

    while (!ifs.eof()) {
        kit::io::load_xyz(ifs, comment, elements, xyz);
        if (elements.empty())
            break;

        iss.clear();
        iss.str(comment);

        double binding_energy;
        double nb_energy;

        iss >> binding_energy >> std::ws
            >> nb_energy ;

        if (iss.fail() and is_training_set) {
            std::ostringstream oss;
            oss << "'" << filename << "' : configuration #"
                << (n_nb_sys + 1) << " : unexpected text '"
                << iss.str() << "' instead of total/twobody/onebody[2] energies";
            throw std::runtime_error(oss.str());
        } else {
            binding_energy = 0.0;
            nb_energy = 0.0;
        }

        ts.push_back(nb_system());

        ts.rbegin()->binding_energy = binding_energy;
        ts.rbegin()->nb_energy = nb_energy;

        ts.rbegin()->xyz = xyz;

        ++n_nb_sys;
    }

    return n_nb_sys;
}

} // namespace tset
