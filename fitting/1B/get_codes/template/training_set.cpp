#include <cassert>

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "io-xyz.h"
#include "training_set.h"
//#include "xyz-water-utils.h"

namespace tset {

size_t load_monomers(const char* filename, std::vector<monomer>& ts)
{
    assert(filename);

    std::ifstream ifs(filename);
    if (!ifs) {
        std::ostringstream oss;
        oss << "could not open '" << filename << "' for reading";
        throw std::runtime_error(oss.str());
    }

    size_t nmonomers(0);

    std::string comment;
    std::vector<std::string> elements;
    std::vector<double> xyz;

    std::istringstream iss;

    while (!ifs.eof()) {
        kit::io::load_xyz(ifs, comment, elements, xyz);
        if (elements.empty())
            break;
/*
        if (elements.size() != 9 || !kit::is_water(elements, xyz, false)) {
            std::ostringstream oss;
            oss << "'" << filename << "' : configuration #"
                << (nmonomers + 1) << " is not a water monomer";
            throw std::runtime_error(oss.str());
        }
*/
        iss.clear();
        iss.str(comment);

        double energy_onebody;

        iss >> energy_onebody; 

        if (iss.fail()) {
            std::ostringstream oss;
            oss << "'" << filename << "' : configuration #"
                << (nmonomers + 1) << " : unexpected text '"
                << iss.str() << "' instead of the three-body energy  ";
            throw std::runtime_error(oss.str());
        }

        ts.push_back(monomer());

        ts.rbegin()->energy_onebody = energy_onebody;

        std::copy(xyz.begin(), xyz.end(), ts.rbegin()->xyz);

        ++nmonomers;
    }

    return nmonomers;
}

} // namespace tset
