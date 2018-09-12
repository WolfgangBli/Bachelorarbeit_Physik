#ifndef COMPLEXVECTOR
#define COMPLEXVECTOR

#include <cmath>
#include <complex>
#include <vector>

#endif

#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <boost/math/special_functions/bessel.hpp> //Besselfunktion

std::vector<std::complex<double> > chebyshev(std::vector<std::complex<double> > psi,double t,int L, int N);

#endif
