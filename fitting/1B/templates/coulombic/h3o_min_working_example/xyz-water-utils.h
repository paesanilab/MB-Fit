#ifndef XYZ_WATER_UTILS_H
#define XYZ_WATER_UTILS_H

#include <string>
#include <vector>

namespace kit {

bool is_water(std::vector<std::string>& elements,
              std::vector<double>&      xyz,
              bool reorder = true, const double& rOH = 1.3);

void make_water_elements(size_t nw, std::vector<std::string>& elements);

} // namespace kit

#endif // XYZ_WATER_UTILS_H
