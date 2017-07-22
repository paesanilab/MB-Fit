
#ifndef DISPERSION_H
#define DISPERSION_H

#include "tang-toennies.h"
#include <cmath>
#include <algorithm>



struct x2b_disp {
  x2b_disp();
  x2b_disp(double *, double * , size_t, size_t);
  ~x2b_disp();

  const double m_C6_AA = 310.915 ; 
  const double m_C6_AB = 216.702 ; 
  const double m_C6_BB = 172.445 ; 
  const double m_d6_AA = 3.20676 ; 
  const double m_d6_AB = 3.38985 ; 
  const double m_d6_BB = 3.74833 ; 


  const double m_C8 = 0.0;
  const double m_d8 = 0.0;

  const double if6 = 1.0/x2o::factorial<6>();
  const double if8 = 1.0/x2o::factorial<8>();

  double * xyz1;
  double * xyz2;

  double get_dispersion();
  double get_dispersion(double * grdx);


  inline double x6(const double& C6, const double& d6,
                  const double& C8, const double& d8,
                  const double* p1, const double* p2,
                        double* g1,       double* g2)
  {
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double d6r = d6*r;
    const double tt6 = x2o::tang_toennies(6, d6r);

    const double d8r = d8*r;
    const double tt8 = x2o::tang_toennies(8, d8r);


    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;
    const double inv_r8 = inv_r6*inv_rsq;

    const double e6 = C6*tt6*inv_r6;
    const double e8 = C8*tt8*inv_r8;

    const double grd = (6*e6 + 8*e8)*inv_rsq
        - (C6*std::pow(d6, 7)*if6*std::exp(-d6r)
        +  C8*std::pow(d8, 9)*if8*std::exp(-d8r))/r;

    g1[0] += dx*grd;
    g2[0] -= dx*grd;

    g1[1] += dy*grd;
    g2[1] -= dy*grd;

    g1[2] += dz*grd;
    g2[2] -= dz*grd;

    return - (e6 + e8);
  }

  inline double x6(const double& C6, const double& d6,
                  const double& C8, const double& d8,
                  const double* p1, const double* p2)
  {
    const double dx = p1[0] - p2[0];
    const double dy = p1[1] - p2[1];
    const double dz = p1[2] - p2[2];

    const double rsq = dx*dx + dy*dy + dz*dz;
    const double r = std::sqrt(rsq);

    const double d6r = d6*r;
    const double tt6 = x2o::tang_toennies(6, d6r);

    const double d8r = d8*r;
    const double tt8 = x2o::tang_toennies(8, d8r);


    const double inv_rsq = 1.0/rsq;
    const double inv_r6 = inv_rsq*inv_rsq*inv_rsq;
    const double inv_r8 = inv_r6*inv_rsq;

    const double e6 = C6*tt6*inv_r6;
    const double e8 = C8*tt8*inv_r8;

    return - (e6 + e8);
  }
  
};

#endif

