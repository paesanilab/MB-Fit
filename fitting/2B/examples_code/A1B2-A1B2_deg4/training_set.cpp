#include <cassert>

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "io-xyz.h"
#include "training_set.h"

namespace tset {

size_t load_dimers(const char* filename, std::vector<dimer>& ts)
{
    assert(filename);

    std::ifstream ifs(filename);
    if (!ifs) {
        std::ostringstream oss;
        oss << "could not open '" << filename << "' for reading";
        throw std::runtime_error(oss.str());
    }

    size_t ndimers(0);

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

        double energy_total;
        double energy_twobody;
        double energy_onebody[2];

        iss >> energy_total >> std::ws
            >> energy_twobody >> std::ws
            >> energy_onebody[0] >> std::ws
            >> energy_onebody[1];

        if (iss.fail()) {
            std::ostringstream oss;
            oss << "'" << filename << "' : configuration #"
                << (ndimers + 1) << " : unexpected text '"
                << iss << "' instead of total/twobody/onebody[2] energies";
            throw std::runtime_error(oss.str());
        }

        ts.push_back(dimer());

        ts.rbegin()->energy_total = energy_total;
        ts.rbegin()->energy_twobody = energy_twobody;
        ts.rbegin()->energy_onebody[0] = energy_onebody[0];
        ts.rbegin()->energy_onebody[1] = energy_onebody[1];

        std::copy(xyz.begin(), xyz.end(), ts.rbegin()->xyz);

        ++ndimers;
    }

    return ndimers;
}

} // namespace tset
