#ifndef VSITES_H
#define VSITES_H

struct vsites {
    //void TwoParticleAverageSite() {}
    //void ThreeParticleAverageSite() {}
    void OutOfPlaneSite(const double& w12, const double& w13,
                        const double& wcross, const double x1[3],
                        const double y1[3], const double y2[3],
                        double vs[3]);
    //void LocalCoordinatesSite{}
};

#endif // VSITES_H
