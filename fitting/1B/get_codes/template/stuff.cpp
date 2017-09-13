#include <cmath>
#include <cassert>

#include "stuff.h"

namespace x2o {

void align_water_dimer(double* xyz)
{
    assert(xyz);

    // put first atom at the origin

    {
        const double orig[3] = {xyz[0], xyz[1], xyz[2]};
        for (unsigned n = 0; n < 18; ++n)
            xyz[n] -= orig[n%3];
    }

    // rotate so that the 1st monomer is in the XY plane with OH along X

    double ex[3] = {xyz[3], xyz[4], xyz[5]};
    x2o::normalize(ex);

    double ey[3] = {0.5*(xyz[3] + xyz[6]),
                    0.5*(xyz[4] + xyz[7]),
                    0.5*(xyz[5] + xyz[8])};

    {
        const double tmp = ex[0]*ey[0] + ex[1]*ey[1] + ex[2]*ey[2];

        ey[0] -= tmp*ex[0];
        ey[1] -= tmp*ex[1];
        ey[2] -= tmp*ex[2];
    }

    x2o::normalize(ey);

    double ez[3] = {ex[1]*ey[2] - ex[2]*ey[1],
                    ex[2]*ey[0] - ex[0]*ey[2],
                    ex[0]*ey[1] - ex[1]*ey[0]};

    const double R[3][3] = {
        {ex[0], ex[1], ex[2]},
        {ey[0], ey[1], ey[2]},
        {ez[0], ez[1], ez[2]}
    };

    for (unsigned n = 0; n < 18; n += 3) {
        const double orig[3] = {xyz[n], xyz[n + 1], xyz[n + 2]};
        xyz[n] = xyz[n + 1] = xyz[n + 2] = 0.0;
        for (unsigned i = 0; i < 3; ++i)
            for (unsigned j = 0; j < 3; ++j)
                xyz[n + i] += R[i][j]*orig[j];
    }
}

} // namespace x2o
