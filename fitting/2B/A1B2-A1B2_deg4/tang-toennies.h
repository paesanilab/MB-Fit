#ifndef TANG_TOENNIES_H
#define TANG_TOENNIES_H

#include <cstddef>

namespace x2o {

//
//  T(n, x) = 1 - exp(-x)*(1 + x/1! + x^2/2! + x^3/3! + ... + x^n/n!)
//  diff(T(n, x), x) = T(n - 1, x) - T(n, x) = exp(-x)*x^n/n!
//

double tang_toennies(int, const double&);

template <int N>
inline int factorial()
{
    return N*factorial<N-1>();
}

template<>
inline int factorial<0>()
{
    return 1;
}

} // namespace x2o

#endif // TANG_TOENNIES_H
