#ifndef STUFF_H
#define STUFF_H

#include <cmath>
#include <cstddef>

////////////////////////////////////////////////////////////////////////////////

namespace x2o {

////////////////////////////////////////////////////////////////////////////////

inline double distance(const double* p1, const double* p2)
{
    double result(0);

    for (size_t k = 0; k < 3; ++k) {
        const double delta = p1[k] - p2[k];
        result += delta*delta;
    }

    return std::sqrt(result);
}

//----------------------------------------------------------------------------//

inline double angle(const double* p1, const double* p2, const double* p3)
{
    const double v1[3] = {p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]};
    const double v2[3] = {p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]};

    const double v1v2 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    const double v1sq = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
    const double v2sq = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];

    return 180.0*std::acos(v1v2/std::sqrt(v1sq*v2sq))/M_PI;
}

//----------------------------------------------------------------------------//

inline void normalize(double* v)
{
    const double vv = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

    v[0] /= vv;
    v[1] /= vv;
    v[2] /= vv;
}

//----------------------------------------------------------------------------//

void align_water_dimer(double*);

////////////////////////////////////////////////////////////////////////////////

} // namespace x2o

////////////////////////////////////////////////////////////////////////////////

#endif // STUFF_H
