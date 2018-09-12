#ifndef COMPLEXVECTOR
#define COMPLEXVECTOR

#include <cmath>
#include <complex>
#include <vector>

#endif

#ifndef KRYLOV_H
#define KRYLOV_H

#include <lapacke.h>
#include <cblas.h>

std::vector<std::complex<double> > krylov(std::vector<std::complex<double> > psi,double t,int L, int m);

#endif
