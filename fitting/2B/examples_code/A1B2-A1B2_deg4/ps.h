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

#ifndef H2O_PS_H
#define H2O_PS_H

namespace h2o { namespace ps {

double pot_nasa(const double*, double*);
void dms_nasa(const double&, const double&, const double&,
              const double*, double* q3, double* dq3, bool ttm3);

}} // namespace h2o::ps

#endif // H2O_PS_H
