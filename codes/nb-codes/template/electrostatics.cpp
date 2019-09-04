// Copyright 2013 Volodymyr Babin <vb27606@gmail.com>
//
// This is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// The code is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You can find a copy of the GNU General Public License at
// http://www.gnu.org/licenses/.

#include <cmath>
#include <cassert>

#include <algorithm>
#include <stdexcept>

#include "constants.h"
#include "electrostatics.h"

//#define DEBUG yes

#ifdef DEBUG
#  include <iomanip>
#  include <iostream>
#endif // DEBUG

#define CHOLESKY yes

const double DEBYE = constants::DEBYE;
const double CHARGECON = constants::CHARGECON;
namespace ttm {

electrostatics::electrostatics() //constructor with no arguments 
: m_size(0), m_converged(false)
{
    allocate(0);
}

electrostatics::~electrostatics() //destructor
{
    allocate(0);
}

void electrostatics::allocate(size_t n) 
{
    if (m_size == n) {
        return;
	}

    if (m_size > 0) {
		delete[] m_memory;
	}

    m_size = n;

    if (n > 0)
        m_memory = new double[n*(6*3 + 1 + 9*n)];
}

size_t electrostatics::size() const
{
    return m_size;
}

const double* electrostatics::dipoles() const

{
    return m_memory;
}

bool electrostatics::dipoles_converged() const
{
    return m_converged;
}

const double* electrostatics::phi() const
{
    return m_memory + 3*m_size + 3*m_size; // see operator() below
}

double electrostatics::operator()(size_t natom, const double* charge,
    const double* polfac, const double* polar, const double* xyz,
    const electrostatics::excluded_set_type& excluded12,
    const electrostatics::excluded_set_type& excluded13,
    const electrostatics::excluded_set_type& excluded14,
    const int * is_w,
    const smear& smr, double* grad)
{
    allocate(natom);

    double* dipole      = m_memory;

#   ifdef CHOLESKY  
    double* polar_sqrt  = dipole + 3*natom;
    double* phi         = polar_sqrt + 3*natom;
#   else
    double* dipole_prev = dipole + 3*natom;
    double* phi         = dipole_prev + 3*natom;
#   endif
    double* Efq         = phi + natom;
    double* Efd         = Efq + 3*natom;
    double* ddt         = Efd + 3*natom;

    std::fill(phi, phi + natom*(1 + 2*3 + 3*3*natom), 0.0);

#   ifdef CHOLESKY
    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        const double sqrtalpha = std::sqrt(polar[i]);

        for (size_t k = 0; k < 3; ++k)
            polar_sqrt[i3 + k] = sqrtalpha;
    }
#   endif

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;

        for (size_t j = i + 1; j < natom; ++j) {
            const size_t j3 = 3*j;

            double Rij[3], Rsq(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
            }

            // charge-charge

            const bool ij_12 =
                (excluded12.find(std::make_pair(i, j)) != excluded12.end());

            const bool ij_13 =
                (excluded13.find(std::make_pair(i, j)) != excluded13.end());

            const bool ij_14 =
                (excluded14.find(std::make_pair(i, j)) != excluded14.end());


            const bool ij_bonded = ij_12 || ij_13 || ij_14;

            if (!ij_bonded) {
                double ts0, ts1;
                smr.smear01(Rsq, polfac[i]*polfac[j], smr.aCC(), ts0, ts1);

                phi[i] += ts0*charge[j];
                phi[j] += ts0*charge[i];
                for (size_t k = 0; k < 3; ++k) {
                    Efq[i3 + k] += ts1*charge[j]*Rij[k];
                    Efq[j3 + k] -= ts1*charge[i]*Rij[k];
                }
            }

	    // dipole-dipole tensor

            double aDD = ij_bonded ?
             (ij_12 ? smr.aDD_intra_12() : ij_13 ? smr.aDD_intra_13()
                                         : smr.aDD_intra_14())
                                         : smr.aDD_inter();
            if (ij_bonded && is_w[i] == 1) {
              if (ij_12) {
                aDD = 0.626;
              } else {
                aDD = 0.055;
              }
            }
//DEBUG
//	    std::cout << " i = " << i << "    j = " << j << 
//		"   aDD = " << aDD << std::endl;

            double ts1, ts2;
            smr.smear2(Rsq, polfac[i]*polfac[j], aDD, ts1, ts2);

            double dd3[3][3];

            dd3[0][0] = 3.0*ts2*Rij[0]*Rij[0] - ts1;
            dd3[1][1] = 3.0*ts2*Rij[1]*Rij[1] - ts1;
            dd3[2][2] = 3.0*ts2*Rij[2]*Rij[2] - ts1;
            dd3[0][1] = 3.0*ts2*Rij[0]*Rij[1];
            dd3[0][2] = 3.0*ts2*Rij[0]*Rij[2];
            dd3[1][2] = 3.0*ts2*Rij[1]*Rij[2];
            dd3[1][0] = dd3[0][1];
            dd3[2][0] = dd3[0][2];
            dd3[2][1] = dd3[1][2];

#           ifdef CHOLESKY
            const double aiaj = polar_sqrt[i3]*polar_sqrt[j3];
            for (size_t k = 0; k < 3; ++k)
                for (size_t l = 0; l < 3; ++l)
                    ddt[3*natom*(i3 + k) + j3 + l] = -aiaj*dd3[k][l];
#           else
            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    ddt[3*natom*(i3 + k) + j3 + l] = dd3[k][l];
                    ddt[3*natom*(j3 + k) + i3 + l] = dd3[k][l];
                }
            }
#           endif // CHOLESKY
        }
    }


#   ifdef CHOLESKY

    // populate the diagonal
    for (size_t i = 0; i < 3*natom; ++i)
        ddt[3*natom*i + i] = 1.0;

    // perform Cholesky decomposition ddt = L*L^T
    // storing L in the lower triangle of ddt and
    // its diagonal in Efd

    const size_t natom3 = 3*natom;

    for (size_t i = 0; i < natom3; ++i) {
        for (size_t j = i; j < natom3; ++j) {

            double sum = ddt[i*natom3 + j];
            for (size_t k = 0; k < i; ++k)
                sum -= ddt[i*natom3 + k]*ddt[j*natom3 + k];

            if(i == j) {
                if (!(sum > 0.0))
                    throw std::logic_error("induced dipoles diverged");
                Efd[i] = std::sqrt(sum);
            } else {
                ddt[j*natom3 + i] = sum/Efd[i];
            }
        }
    }

    // solve L*y = sqrt(a)*Efq storing y in dipole[]
    for (size_t i = 0; i < natom3; ++i) {
        double sum = polar_sqrt[i]*Efq[i];
        for (size_t k = 0; k < i; ++k)
            sum -= ddt[i*natom3 + k]*dipole[k];

        dipole[i] = sum/Efd[i];
    }

    // solve L^T*x = y
    for (size_t i = natom3; i > 0; --i) {
        const size_t i1 = i - 1;
        double sum = dipole[i1];
        for (size_t k = i; k < natom3; ++k)
            sum -= ddt[k*natom3 + i1]*dipole[k];

        dipole[i1] = sum/Efd[i1];
    }

    for (size_t i = 0; i < natom3; ++i)
        dipole[i] *= polar_sqrt[i];

    m_converged = true;

#   else

    // Calculate induced dipoles iteratively

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        for (size_t k = 0; k < 3; ++k) {
            dipole[i3 + k] = polar[k]*Efq[i3 + k];
	    // FIXME: k should be i!?
            //dipole[i3 + k] = polar[i]*Efq[i3 + k];
            dipole_prev[i3 + k] = dipole[i3 + k];
        }
    }

    const double dmix = 0.7;
    const double stath = DEBYE/CHARGECON/(std::sqrt(double(natom)));

    m_converged = false;

    for (size_t iter = 0; iter < dip_maxiter; ++iter) {

        for (size_t k = 0; k < 3*natom; ++k)
            for (size_t l = 0; l < 3*natom; ++l)
                Efd[k] += ddt[3*natom*k + l]*dipole[l];

        double deltadip = 0.0;

        for (size_t i = 0; i < natom; ++i) {
            const size_t i3 = 3*i;
            for(size_t k = 0; k < 3; ++k) {
                dipole[i3 + k] = polar[i]*(Efq[i3 + k] + Efd[i3 + k]);
                dipole[i3 + k] = dmix*dipole[i3 + k]
                    + (1.0 - dmix)*dipole_prev[i3 + k];

                const double delta = dipole[i3 + k] - dipole_prev[i3 + k];
                deltadip += delta*delta;
            }
        }

        deltadip = std::sqrt(deltadip)*stath;

        if (deltadip < dip_tolerance) {
            m_converged = true;
            break; // converged!
        } else {
            std::copy(dipole, dipole + 3*natom, dipole_prev);
            std::fill(Efd, Efd + 3*natom, 0.0);
        }
    }

#   endif // CHOLESKY

    double Eelec = 0.0;
    for (size_t i = 0; i < natom; ++i)
        Eelec += phi[i]*charge[i];

    Eelec *= 0.5;

    double Eind = 0.0;
    for (size_t i = 0; i < 3*natom; ++i)
        Eind -= dipole[i]*Efq[i];

    Eind *= 0.5;

    m_Eind = Eind;
    m_Eelec = Eelec;

#   ifdef DEBUG
    double U_spring(0);
    for (size_t n = 0; n < natom; ++n) {
        if (polar[n] > 0.0)
            for (size_t k = 0; k < 3; ++k)
                U_spring += dipole[3*n + k]*dipole[3*n + k]/polar[n];
    }

    U_spring /= 2;
#   endif // DEBUG

    // that is all if no gradients requested
    if (grad == 0) {
        return Eelec + Eind;
    }

    // charge-charge interactions
    for (size_t i = 0; i < natom; ++i)
        for (size_t k = 0; k < 3; ++k)
            grad[3*i + k] = -charge[i]*Efq[3*i + k];

#   ifdef DEBUG
    double U_CD(0);
    double U_DD(0);
#   endif // DEBUG

    // charge-dipole interactions
    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        for (size_t j = 0; j < natom; ++j) {
            const size_t j3 = 3*j;
            const double qj = charge[j];

            const bool skip_ij = (i == j)
            || (i < j ? (excluded12.find(std::make_pair(i, j)) != excluded12.end()
                      || excluded13.find(std::make_pair(i, j)) != excluded13.end()
                      || excluded14.find(std::make_pair(i, j)) != excluded14.end())
                      : (excluded12.find(std::make_pair(j, i)) != excluded12.end()
                      || excluded13.find(std::make_pair(j, i)) != excluded13.end()
                      || excluded14.find(std::make_pair(j, i)) != excluded14.end()));

            if (skip_ij)
                continue; // skip this (i, j) pair

            double Rij[3], Rsq(0), diR(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
                diR += dipole[i3 + k]*Rij[k];
            }

            double ts1, ts2;
            smr.smear2(Rsq, polfac[i]*polfac[j], smr.aCD(), ts1, ts2);

            for (size_t k = 0; k < 3; ++k) {
                const double derij =
                    qj*(3*ts2*diR*Rij[k] - ts1*dipole[i3 + k]);

                grad[i3 + k] += derij;
                grad[j3 + k] -= derij;
            }

            phi[j] -= ts1*diR;

#           ifdef DEBUG
            U_CD -= qj*diR*ts1; // sign is due to Rij
#           endif // DEBUG

        }
    }

    // dipole-dipole interactions
    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        for (size_t j = i + 1; j < natom; ++j) {
            const size_t j3 = 3*j;

            double Rij[3], Rsq(0), diR(0), djR(0), didj(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
                diR += dipole[i3 + k]*Rij[k];
                djR += dipole[j3 + k]*Rij[k];
                didj += dipole[i3 + k]*dipole[j3 + k];
            }

            const bool ij_12 =
                excluded12.find(std::make_pair(i, j)) != excluded12.end();

            const bool ij_13 =
                excluded13.find(std::make_pair(i, j)) != excluded13.end();

            const bool ij_14 =
                excluded14.find(std::make_pair(i, j)) != excluded14.end();

            const bool ij_bonded = ij_12 || ij_13 || ij_14;

            double aDD = ij_bonded ?
             (ij_12 ? smr.aDD_intra_12() : ij_13 ? smr.aDD_intra_13()
                                         : smr.aDD_intra_14())
                                         : smr.aDD_inter();
            if (ij_bonded && is_w[i] == 1) {
              if (ij_12) {
                aDD = 0.626;
              } else {
                aDD = 0.055;
              }
            }

            double ts1, ts2, ts3;
            smr.smear3(Rsq, polfac[i]*polfac[j], aDD, ts1, ts2, ts3);

            for (size_t k = 0; k < 3; ++k) {
                const double derij =
                   - 3*ts2*(didj*Rij[k] + djR*dipole[i3 + k]
                                        + diR*dipole[j3 + k])
                   + 15*ts3*diR*djR*Rij[k];

                grad[i3 + k] += derij;
                grad[j3 + k] -= derij;
            }
#           ifdef DEBUG
            U_DD += didj*ts1 - 3*diR*djR*ts2;
#           endif // DEBUG
        }
    }

#   ifdef DEBUG
    std::cerr << "    U_CC = " << Eelec << std::endl;
    std::cerr << "    U_CD = " << U_CD << std::endl;
    std::cerr << "    U_DD = " << U_DD << std::endl;
    std::cerr << "U_spring = " << U_spring << std::endl;

    std::cerr << "\nU_CD + U_DD + U_spring = "
              << (U_CD + U_DD + U_spring) << std::endl;
    std::cerr << "E_ind = " << Eind << std::endl;
#   endif // DEBUG

    return Eelec + Eind;
}

//The other operator which returns Eelec and Eind separately
double electrostatics::operator()(size_t natom, const double* charge,
    const double* polfac, const double* polar, const double* xyz,
    const electrostatics::excluded_set_type& excluded12,
    const electrostatics::excluded_set_type& excluded13,
    const electrostatics::excluded_set_type& excluded14,
    const int * is_w,
    const smear& smr, double& E_elec, double& E_ind, double* grad)
{
    allocate(natom);

    double* dipole      = m_memory;

#   ifdef CHOLESKY  
    double* polar_sqrt  = dipole + 3*natom;
    double* phi         = polar_sqrt + 3*natom;
#   else
    double* dipole_prev = dipole + 3*natom;
    double* phi         = dipole_prev + 3*natom;
#   endif
    double* Efq         = phi + natom;
    double* Efd         = Efq + 3*natom;
    double* ddt         = Efd + 3*natom;

    std::fill(phi, phi + natom*(1 + 2*3 + 3*3*natom), 0.0);

#   ifdef CHOLESKY
    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        const double sqrtalpha = std::sqrt(polar[i]);

        for (size_t k = 0; k < 3; ++k)
            polar_sqrt[i3 + k] = sqrtalpha;
    }
#   endif

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;

        for (size_t j = i + 1; j < natom; ++j) {
            const size_t j3 = 3*j;

            double Rij[3], Rsq(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
            }

            // charge-charge

            const bool ij_12 =
                (excluded12.find(std::make_pair(i, j)) != excluded12.end());

            const bool ij_13 =
                (excluded13.find(std::make_pair(i, j)) != excluded13.end());

            const bool ij_14 =
                (excluded14.find(std::make_pair(i, j)) != excluded14.end());

            const bool ij_bonded = ij_12 || ij_13 || ij_14;

            if (!ij_bonded) {
                double ts0, ts1;
                smr.smear01(Rsq, polfac[i]*polfac[j], smr.aCC(), ts0, ts1);

                phi[i] += ts0*charge[j];
                phi[j] += ts0*charge[i];
                for (size_t k = 0; k < 3; ++k) {
                    Efq[i3 + k] += ts1*charge[j]*Rij[k];
                    Efq[j3 + k] -= ts1*charge[i]*Rij[k];
                }
            }

	    // dipole-dipole tensor

            double aDD = ij_bonded ?
             (ij_12 ? smr.aDD_intra_12() : ij_13 ? smr.aDD_intra_13()
                                         : smr.aDD_intra_14())
                                         : smr.aDD_inter();
            if (ij_bonded && is_w[i] == 1) {
              if (ij_12) {
                aDD = 0.626;
              } else {
                aDD = 0.055;
              }
            }

            double ts1, ts2;
            smr.smear2(Rsq, polfac[i]*polfac[j], aDD, ts1, ts2);

            double dd3[3][3];

            dd3[0][0] = 3.0*ts2*Rij[0]*Rij[0] - ts1;
            dd3[1][1] = 3.0*ts2*Rij[1]*Rij[1] - ts1;
            dd3[2][2] = 3.0*ts2*Rij[2]*Rij[2] - ts1;
            dd3[0][1] = 3.0*ts2*Rij[0]*Rij[1];
            dd3[0][2] = 3.0*ts2*Rij[0]*Rij[2];
            dd3[1][2] = 3.0*ts2*Rij[1]*Rij[2];
            dd3[1][0] = dd3[0][1];
            dd3[2][0] = dd3[0][2];
            dd3[2][1] = dd3[1][2];

#           ifdef CHOLESKY
            const double aiaj = polar_sqrt[i3]*polar_sqrt[j3];
            for (size_t k = 0; k < 3; ++k)
                for (size_t l = 0; l < 3; ++l)
                    ddt[3*natom*(i3 + k) + j3 + l] = -aiaj*dd3[k][l];
#           else
            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    ddt[3*natom*(i3 + k) + j3 + l] = dd3[k][l];
                    ddt[3*natom*(j3 + k) + i3 + l] = dd3[k][l];
                }
            }
#           endif // CHOLESKY
        }
    }


#   ifdef CHOLESKY

    // populate the diagonal
    for (size_t i = 0; i < 3*natom; ++i)
        ddt[3*natom*i + i] = 1.0;

    // perform Cholesky decomposition ddt = L*L^T
    // storing L in the lower triangle of ddt and
    // its diagonal in Efd

    const size_t natom3 = 3*natom;

    for (size_t i = 0; i < natom3; ++i) {
        for (size_t j = i; j < natom3; ++j) {

            double sum = ddt[i*natom3 + j];
            for (size_t k = 0; k < i; ++k)
                sum -= ddt[i*natom3 + k]*ddt[j*natom3 + k];

            if(i == j) {
                if (!(sum > 0.0))
                    throw std::logic_error("induced dipoles diverged");
                Efd[i] = std::sqrt(sum);
            } else {
                ddt[j*natom3 + i] = sum/Efd[i];
            }
        }
    }

    // solve L*y = sqrt(a)*Efq storing y in dipole[]
    for (size_t i = 0; i < natom3; ++i) {
        double sum = polar_sqrt[i]*Efq[i];
        for (size_t k = 0; k < i; ++k)
            sum -= ddt[i*natom3 + k]*dipole[k];

        dipole[i] = sum/Efd[i];
    }

    // solve L^T*x = y
    for (size_t i = natom3; i > 0; --i) {
        const size_t i1 = i - 1;
        double sum = dipole[i1];
        for (size_t k = i; k < natom3; ++k)
            sum -= ddt[k*natom3 + i1]*dipole[k];

        dipole[i1] = sum/Efd[i1];
    }

    for (size_t i = 0; i < natom3; ++i)
        dipole[i] *= polar_sqrt[i];

    m_converged = true;

#   else

    // Calculate induced dipoles iteratively

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        for (size_t k = 0; k < 3; ++k) {
            dipole[i3 + k] = polar[k]*Efq[i3 + k];
	    // FIXME: k should be i!?
            //dipole[i3 + k] = polar[i]*Efq[i3 + k];
            dipole_prev[i3 + k] = dipole[i3 + k];
        }
    }

    const double dmix = 0.7;
    const double stath = DEBYE/CHARGECON/(std::sqrt(double(natom)));

    m_converged = false;

    for (size_t iter = 0; iter < dip_maxiter; ++iter) {

        for (size_t k = 0; k < 3*natom; ++k)
            for (size_t l = 0; l < 3*natom; ++l)
                Efd[k] += ddt[3*natom*k + l]*dipole[l];

        double deltadip = 0.0;

        for (size_t i = 0; i < natom; ++i) {
            const size_t i3 = 3*i;
            for(size_t k = 0; k < 3; ++k) {
                dipole[i3 + k] = polar[i]*(Efq[i3 + k] + Efd[i3 + k]);
                dipole[i3 + k] = dmix*dipole[i3 + k]
                    + (1.0 - dmix)*dipole_prev[i3 + k];

                const double delta = dipole[i3 + k] - dipole_prev[i3 + k];
                deltadip += delta*delta;
            }
        }

        deltadip = std::sqrt(deltadip)*stath;

        if (deltadip < dip_tolerance) {
            m_converged = true;
            break; // converged!
        } else {
            std::copy(dipole, dipole + 3*natom, dipole_prev);
            std::fill(Efd, Efd + 3*natom, 0.0);
        }
    }

#   endif // CHOLESKY

    double Eelec = 0.0;
    for (size_t i = 0; i < natom; ++i)
        Eelec += phi[i]*charge[i];

    Eelec *= 0.5;

    double Eind = 0.0;
    for (size_t i = 0; i < 3*natom; ++i)
        Eind -= dipole[i]*Efq[i];

    Eind *= 0.5;

    m_Eind = Eind;
    m_Eelec = Eelec;

#   ifdef DEBUG
    double U_spring(0);
    for (size_t n = 0; n < natom; ++n) {
        if (polar[n] > 0.0)
            for (size_t k = 0; k < 3; ++k)
                U_spring += dipole[3*n + k]*dipole[3*n + k]/polar[n];
    }

    U_spring /= 2;
#   endif // DEBUG

    // that is all if no gradients requested
    if (grad == 0) {
	E_elec = Eelec;
	E_ind = Eind;
	return Eelec + Eind;
    }

    // charge-charge interactions
    for (size_t i = 0; i < natom; ++i)
        for (size_t k = 0; k < 3; ++k)
            grad[3*i + k] = -charge[i]*Efq[3*i + k];

#   ifdef DEBUG
    double U_CD(0);
    double U_DD(0);
#   endif // DEBUG

    // charge-dipole interactions
    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        for (size_t j = 0; j < natom; ++j) {
            const size_t j3 = 3*j;
            const double qj = charge[j];

            const bool skip_ij = (i == j)
            || (i < j ? (excluded12.find(std::make_pair(i, j)) != excluded12.end()
                      || excluded13.find(std::make_pair(i, j)) != excluded13.end()
                      || excluded14.find(std::make_pair(i, j)) != excluded14.end())
                      : (excluded12.find(std::make_pair(j, i)) != excluded12.end()
                      || excluded13.find(std::make_pair(j, i)) != excluded13.end()
                      || excluded14.find(std::make_pair(j, i)) != excluded14.end()));

            if (skip_ij)
                continue; // skip this (i, j) pair

            double Rij[3], Rsq(0), diR(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
                diR += dipole[i3 + k]*Rij[k];
            }

            double ts1, ts2;
            smr.smear2(Rsq, polfac[i]*polfac[j], smr.aCD(), ts1, ts2);

            for (size_t k = 0; k < 3; ++k) {
                const double derij =
                    qj*(3*ts2*diR*Rij[k] - ts1*dipole[i3 + k]);

                grad[i3 + k] += derij;
                grad[j3 + k] -= derij;
            }

            phi[j] -= ts1*diR;

#           ifdef DEBUG
            U_CD -= qj*diR*ts1; // sign is due to Rij
#           endif // DEBUG

        }
    }

    // dipole-dipole interactions
    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        for (size_t j = i + 1; j < natom; ++j) {
            const size_t j3 = 3*j;

            double Rij[3], Rsq(0), diR(0), djR(0), didj(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
                diR += dipole[i3 + k]*Rij[k];
                djR += dipole[j3 + k]*Rij[k];
                didj += dipole[i3 + k]*dipole[j3 + k];
            }

            const bool ij_12 =
                excluded12.find(std::make_pair(i, j)) != excluded12.end();

            const bool ij_13 =
                excluded13.find(std::make_pair(i, j)) != excluded13.end();

            const bool ij_14 =
                excluded14.find(std::make_pair(i, j)) != excluded14.end();

            const bool ij_bonded = ij_12 || ij_13 || ij_14;

            double aDD = ij_bonded ?
             (ij_12 ? smr.aDD_intra_12() : ij_13 ? smr.aDD_intra_13()
                                         : smr.aDD_intra_14())
                                         : smr.aDD_inter();
            if (ij_bonded && is_w[i] == 1) {
              if (ij_12) {
                aDD = 0.626;
              } else {
                aDD = 0.055;
              }
            }

            double ts1, ts2, ts3;
            smr.smear3(Rsq, polfac[i]*polfac[j], aDD, ts1, ts2, ts3);

            for (size_t k = 0; k < 3; ++k) {
                const double derij =
                   - 3*ts2*(didj*Rij[k] + djR*dipole[i3 + k]
                                        + diR*dipole[j3 + k])
                   + 15*ts3*diR*djR*Rij[k];

                grad[i3 + k] += derij;
                grad[j3 + k] -= derij;
            }
#           ifdef DEBUG
            U_DD += didj*ts1 - 3*diR*djR*ts2;
#           endif // DEBUG
        }
    }

#   ifdef DEBUG
    std::cerr << "    U_CC = " << Eelec << std::endl;
    std::cerr << "    U_CD = " << U_CD << std::endl;
    std::cerr << "    U_DD = " << U_DD << std::endl;
    std::cerr << "U_spring = " << U_spring << std::endl;

    std::cerr << "\nU_CD + U_DD + U_spring = "
              << (U_CD + U_DD + U_spring) << std::endl;
    std::cerr << "E_ind = " << Eind << std::endl;
#   endif // DEBUG

    E_elec = Eelec;
    E_ind = Eind;
    return 0;
}


//
// assumes that () was called
//

const double* electrostatics::polarizability_tensor(size_t natom)
{
#   ifndef CHOLESKY
    assert(false); // works only with Cholesky decomposition
#   endif // CHOLESKY

    double* dipole      = m_memory;
    double* polar_sqrt  = dipole + 3*natom;
    double* phi         = polar_sqrt + 3*natom;
    double* Efq         = phi + natom;
    double* Efd         = Efq + 3*natom;
    double* ddt         = Efd + 3*natom;

    // ddt = L*L^T [done in operator()]
    // L in the lower triangle of ddt and
    // its diagonal in Efd

    // we need B = alpha^{1/2}*ddt^{-1}*alpha^{1/2}

    // compute L^{-1} following NR (ch 2.9)

    const size_t natom3 = 3*natom;

    for (size_t i = 0; i < natom3; ++i) {
        ddt[i*natom3 + i] = 1.0/Efd[i];

        for (size_t j = i + 1; j < natom3; ++j) {

            double sum(0);
            for (size_t k = i; k < j; ++k)
                sum -= ddt[j*natom3 + k]*ddt[k*natom3 + i];

            ddt[j*natom3 + i] = sum/Efd[j];
        }
    }

    std::fill(m_pt, m_pt + 9, 0.0);

    for (size_t i = 0; i < natom; ++i) {
        for (size_t j = 0; j < natom; ++j) {
            for (size_t a = 0; a < 3; ++a) {
                for (size_t b = 0; b < 3; ++b) {

                    double sum(0);

                    const size_t ia = 3*i + a;
                    const size_t jb = 3*j + b;

                    for (size_t k = std::max(ia, jb); k < natom3; ++k)
                        sum += ddt[k*natom3 + ia]*ddt[k*natom3 + jb];

                    m_pt[3*a + b] += sum*polar_sqrt[ia]*polar_sqrt[jb];
                }
            }
        }
    }

    return m_pt;
}

#ifdef ENABLE_FIELD
double electrostatics::energy(size_t natom, const double* charge,
    const double* polfac, const double* polar, const double* xyz,
    const electrostatics::excluded_set_type& excluded12,
    const electrostatics::excluded_set_type& excluded13, //?
    const electrostatics::excluded_set_type& excluded14, //?
    const int * is_w,
    const smear& smr, const double F[3]) //see coulomb.cpp for smear struct
{
    allocate(natom);

    double* dipole      = m_memory;

#   ifdef CHOLESKY  
    double* polar_sqrt  = dipole + 3*natom;
    double* phi         = polar_sqrt + 3*natom;
#   else
    double* dipole_prev = dipole + 3*natom;
    double* phi         = dipole_prev + 3*natom;
#   endif
    double* Efq         = phi + natom;
    double* Efd         = Efq + 3*natom;
    double* ddt         = Efd + 3*natom;

    std::fill(phi, phi + natom*(1 + 2*3 + 3*3*natom), 0.0);

#   ifdef CHOLESKY
    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        const double sqrtalpha = std::sqrt(polar[i]);

        for (size_t k = 0; k < 3; ++k)
            polar_sqrt[i3 + k] = sqrtalpha;
    }
#   endif

    double U_CF(0);

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;

        U_CF -= charge[i]*(F[0]*xyz[i3 + 0]
                         + F[1]*xyz[i3 + 1]
                         + F[2]*xyz[i3 + 2]);

        Efq[i3 + 0] += F[0];
        Efq[i3 + 1] += F[1];
        Efq[i3 + 2] += F[2];

        for (size_t j = i + 1; j < natom; ++j) {
            const size_t j3 = 3*j;

            double Rij[3], Rsq(0);
            for (size_t k = 0; k < 3; ++k) {
                Rij[k] = xyz[i3 + k] - xyz[j3 + k];
                Rsq += Rij[k]*Rij[k];
            }

            // charge-charge

            const bool ij_12 =
                (excluded12.find(std::make_pair(i, j)) != excluded12.end());

            const bool ij_13 =
                (excluded13.find(std::make_pair(i, j)) != excluded13.end());

            const bool ij_14 =
                (excluded14.find(std::make_pair(i, j)) != excluded14.end());

            const bool ij_bonded = ij_12 || ij_13 || ij_14;

            if (!ij_bonded) {
                double ts0, ts1;
                smr.smear01(Rsq, polfac[i]*polfac[j], smr.aCC(), ts0, ts1);

                phi[i] += ts0*charge[j];
                phi[j] += ts0*charge[i];

                for (size_t k = 0; k < 3; ++k) {
                    Efq[i3 + k] += ts1*charge[j]*Rij[k];
                    Efq[j3 + k] -= ts1*charge[i]*Rij[k];
                }
            }

            // dipole-dipole tensor

            double aDD = ij_bonded ?
             (ij_12 ? smr.aDD_intra_12() : ij_13 ? smr.aDD_intra_13()
                                         : smr.aDD_intra_14())
                                         : smr.aDD_inter();
            if (ij_bonded && is_w[i] == 1) {
              if (ij_12) {
                aDD = 0.626;
              } else {
                aDD = 0.055;
              }
            }

            double ts1, ts2;
            smr.smear2(Rsq, polfac[i]*polfac[j], aDD, ts1, ts2);

            double dd3[3][3];

            dd3[0][0] = 3.0*ts2*Rij[0]*Rij[0] - ts1;
            dd3[1][1] = 3.0*ts2*Rij[1]*Rij[1] - ts1;
            dd3[2][2] = 3.0*ts2*Rij[2]*Rij[2] - ts1;
            dd3[0][1] = 3.0*ts2*Rij[0]*Rij[1];
            dd3[0][2] = 3.0*ts2*Rij[0]*Rij[2];
            dd3[1][2] = 3.0*ts2*Rij[1]*Rij[2];
            dd3[1][0] = dd3[0][1];
            dd3[2][0] = dd3[0][2];
            dd3[2][1] = dd3[1][2];

#           ifdef CHOLESKY
            const double aiaj = polar_sqrt[i3]*polar_sqrt[j3];
            for (size_t k = 0; k < 3; ++k)
                for (size_t l = 0; l < 3; ++l)
                    ddt[3*natom*(i3 + k) + j3 + l] = -aiaj*dd3[k][l];
#           else
            for (size_t k = 0; k < 3; ++k) {
                for (size_t l = 0; l < 3; ++l) {
                    ddt[3*natom*(i3 + k) + j3 + l] = dd3[k][l];
                    ddt[3*natom*(j3 + k) + i3 + l] = dd3[k][l];
                }
            }
#           endif // CHOLESKY
        }
    }

#   ifdef CHOLESKY

    // populate the diagonal
    for (size_t i = 0; i < 3*natom; ++i)
        ddt[3*natom*i + i] = 1.0;

    // perform Cholesky decomposition ddt = L*L^T
    // storing L in the lower triangle of ddt and
    // its diagonal in Efd

    const size_t natom3 = 3*natom;

    for (size_t i = 0; i < natom3; ++i) {
        for (size_t j = i; j < natom3; ++j) {

            double sum = ddt[i*natom3 + j];
            for (size_t k = 0; k < i; ++k)
                sum -= ddt[i*natom3 + k]*ddt[j*natom3 + k];

            if(i == j) {
                if (!(sum > 0.0))
                    throw std::logic_error("induced dipoles diverged");
                Efd[i] = std::sqrt(sum);
            } else {
                ddt[j*natom3 + i] = sum/Efd[i];
            }
        }
    }

    // solve L*y = sqrt(a)*Efq storing y in dipole[]
    for (size_t i = 0; i < natom3; ++i) {
        double sum = polar_sqrt[i]*Efq[i];
        for (size_t k = 0; k < i; ++k)
            sum -= ddt[i*natom3 + k]*dipole[k];

        dipole[i] = sum/Efd[i];
    }

    // solve L^T*x = y
    for (size_t i = natom3; i > 0; --i) {
        const size_t i1 = i - 1;
        double sum = dipole[i1];
        for (size_t k = i; k < natom3; ++k)
            sum -= ddt[k*natom3 + i1]*dipole[k];

        dipole[i1] = sum/Efd[i1];
    }

    for (size_t i = 0; i < natom3; ++i)
        dipole[i] *= polar_sqrt[i];

    m_converged = true;

#   else

    // Calculate induced dipoles iteratively

    for (size_t i = 0; i < natom; ++i) {
        const size_t i3 = 3*i;
        for (size_t k = 0; k < 3; ++k) {
            dipole[i3 + k] = polar[k]*Efq[i3 + k];
            dipole_prev[i3 + k] = dipole[i3 + k];
        }
    }

    const double dmix = 0.7;
    const double stath = DEBYE/CHARGECON/(std::sqrt(double(natom)));

    m_converged = false;

    for (size_t iter = 0; iter < dip_maxiter; ++iter) {

        for (size_t k = 0; k < 3*natom; ++k)
            for (size_t l = 0; l < 3*natom; ++l)
                Efd[k] += ddt[3*natom*k + l]*dipole[l];

        double deltadip = 0.0;

        for (size_t i = 0; i < natom; ++i) {
            const size_t i3 = 3*i;
            for(size_t k = 0; k < 3; ++k) {
                dipole[i3 + k] = polar[i]*(Efq[i3 + k] + Efd[i3 + k]);
                dipole[i3 + k] = dmix*dipole[i3 + k]
                    + (1.0 - dmix)*dipole_prev[i3 + k];

                const double delta = dipole[i3 + k] - dipole_prev[i3 + k];
                deltadip += delta*delta;
            }
        }

        deltadip = std::sqrt(deltadip)*stath;

        if (deltadip < dip_tolerance) {
            m_converged = true;
            break; // converged!
        } else {
            std::copy(dipole, dipole + 3*natom, dipole_prev);
            std::fill(Efd, Efd + 3*natom, 0.0);
        }
    }

#   endif // CHOLESKY

    double Eelec = 0.0;
    for (size_t i = 0; i < natom; ++i)
        Eelec += phi[i]*charge[i];

    Eelec *= 0.5;

    double Eind = 0.0;
    for (size_t i = 0; i < 3*natom; ++i)
        Eind -= dipole[i]*Efq[i];

    Eind *= 0.5;

    m_Eind = Eind;
    m_Eelec = Eelec;

#   ifdef DEBUG
    double U_spring(0);
    for (size_t n = 0; n < natom; ++n) {
        if (polar[n] > 0.0)
            for (size_t k = 0; k < 3; ++k)
                U_spring += dipole[3*n + k]*dipole[3*n + k]/polar[n];
    }

    U_spring /= 2;
#   endif // DEBUG

    return Eelec + Eind + U_CF;
}
#endif

} // namespace ttm
