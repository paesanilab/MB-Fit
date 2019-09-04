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

#ifndef ELECTROSTATICS_H
#define ELECTROSTATICS_H

#include <set>
#include <utility>
#include <cstddef>

#include "coulomb.h"

namespace ttm {

struct electrostatics {

    electrostatics();
    ~electrostatics();

    typedef std::set<std::pair<size_t, size_t> > excluded_set_type;

    // energy
    double operator()(size_t N,
                      const double* q,      // size=N, charge
                      const double* polfac, // size=N, damping factor
                      const double* polar,  // size=N, polarizability
                      const double* xyz,    // size=3N, coordinates
                      const excluded_set_type& s12,
                      const excluded_set_type& s13,
                      const excluded_set_type& s14,
                      const int * is_w,
                      const smear&,
                      double* grad = 0);

    double operator()(size_t N,
                      const double* q,      // size=N, charge
                      const double* polfac, // size=N, damping factor
                      const double* polar,  // size=N, polarizability
                      const double* xyz,    // size=3N, coordinates
                      const excluded_set_type& s12,
                      const excluded_set_type& s13,
                      const excluded_set_type& s14,
                      const int * is_w,
                      const smear&,
                      double& E_elec,
		      double& E_ind,
		      double* grad = 0);

/*
#   ifdef ENABLE_FIELD
    double energy(size_t N,
                  const double* q,      // size=N, charge
                  const double* polfac, // size=N, damping factor
                  const double* polar,  // size=N, polarizability
                  const double* xyz,    // size=3N, coordinates
                  const excluded_set_type& s12,
                  const excluded_set_type& s13,
                  const smear&,
                  const double F[3]);   // external field
#   endif
*/
    size_t size() const;

    // size=3*N, induced dipoles
    const double* dipoles() const;
    bool dipoles_converged() const;

    const double& E_ind() const { return m_Eind; }
    const double& E_elec() const { return m_Eelec; }

    const double* polarizability_tensor(size_t N); // 3x3, row major

    // size=N (wrong unless gradients were requested)
    const double* phi() const;

    static const size_t dip_maxiter = 5000;
    static constexpr double dip_tolerance = 1.0e-15;

private:
    size_t  m_size;
    double* m_memory;

    void allocate(size_t);

private:
    bool m_converged;

private:
    double m_Eind, m_Eelec;

private:
    double m_pt[9];
};

} // namespace ttm

#endif // ELECTROSTATICS_H
