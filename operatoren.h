#ifndef COMPLEXVECTOR
#define COMPLEXVECTOR

#include <cmath>
#include <complex>
#include <vector>

#endif

#ifndef OPERATOREN_H
#define OPERATOREN_H

std::vector<std::complex<double> > H(std::vector<std::complex<double> > psi,int L);
std::vector<std::complex<double> > S_z(std::vector<std::complex<double> > psi,int L,int i);
std::vector<std::complex<double> > H_rescaled(std::vector<std::complex<double> > psi,int L,double a, double b);

#endif
