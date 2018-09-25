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

#ifndef COULOMB_H
#define COULOMB_H

namespace ttm {

struct smear {
    virtual ~smear() {};

    virtual double aCC() const = 0;
    virtual double aCD() const = 0;
    virtual double aDD_inter() const = 0;
    virtual double aDD_intra_12() const { return aDD_inter(); }
    virtual double aDD_intra_13() const { return aDD_inter(); }
    virtual double aDD_intra_14() const { return aDD_inter(); }

    virtual void smear01(const double&, const double&, const double&,
                         double&, double&) const = 0;

    virtual void smear2(const double&, const double&, const double&,
                        double&, double&) const = 0;

    virtual void smear3(const double&, const double&, const double&,
                        double&, double&, double&) const = 0;
};

struct smear_ttm4x : public smear {

    smear_ttm4x();

    double aCC() const;
    double aCD() const;
    double aDD_inter() const;
    double aDD_intra_12() const;
    double aDD_intra_13() const;
    double aDD_intra_14() const;

    void smear01(const double&, const double&, const double&,
                 double&, double&) const;

    void smear2(const double&, const double&, const double&,
                double&, double&) const;

    void smear3(const double&, const double&, const double&,
                double&, double&, double&) const;

    double m_aDD_inter;
    double m_aDD_intra_12;
    double m_aDD_intra_13;
    double m_aDD_intra_14;
};

} // namespace ttm

#endif // COULOMB_H
