#include "vsites.h"

void vsites::OutOfPlaneSite(const double& w12,
                            const double& w13,
                            const double& wcross,
                            const double x1[3],
                            const double y1[3],
                            const double y2[3],
                            double vs[3]) {
    double r12[3], r13[3];

    for (int i = 0; i < 3; ++i) {
        r12[i] = y1[i] - x1[i];
        r13[i] = y2[i] - x1[i];
    }
                            
    double rc[3];
    rc[0] = r12[1]*r13[2] - r12[2]*r13[1];
    rc[1] = r12[2]*r13[0] - r12[0]*r13[2];
    rc[2] = r12[0]*r13[1] - r12[1]*r13[0];
    
    vs[0] = x1[0] + w12 * r12[0] + w13 * r13[0] + wcross * rc[0];
    vs[1] = x1[1] + w12 * r12[1] + w13 * r13[1] + wcross * rc[1];
    vs[2] = x1[2] + w12 * r12[2] + w13 * r13[2] + wcross * rc[2];
}
