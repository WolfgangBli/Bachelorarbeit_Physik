#include "chebyshev.h"
#include "operatoren.h"

using namespace std;

vector<complex<double> > chebyshev(vector<complex<double> > psi,double t,int L, int N){
  int l=1<<L;
  double lambda_min=-0.75*(L-1),lambda_max=0.25*(L-1),epsilon=0.025;
  double a=(lambda_max-lambda_min)/(double(2.0)-epsilon),b=(lambda_max+lambda_min)/2.0;
  complex<double> bessel=0,faktor=0,i_n=-1.0i,i_minus=-1.0i;
  vector<complex<double> > t_n_1,t_n_2,t_n,psi_t,Ht_n_1;
  t_n_1.resize(l);
  t_n_2.resize(l);
  t_n.resize(l);
  psi_t.resize(l);
  Ht_n_1.resize(l);
  
  //Werte fuer den nullten und ersten Chebyshevvektor setzen
  t_n_2=psi;
  bessel=boost::math::cyl_bessel_j(0,a*t);
  for(int i=0; i<l; ++i){
    psi_t[i]+=bessel*t_n_2[i];
  }
  t_n_1=H_rescaled(t_n_2,L,a,b);
  bessel=-2.0i*boost::math::cyl_bessel_j(1,a*t);
  for(int i=0; i<l; ++i){
    psi_t[i]+=bessel*t_n_1[i];
  }
  //Rekursive Definition der restlichen Chebyshevvektoren
  for(int n=2; n<N; ++n){
    Ht_n_1=H_rescaled(t_n_1,L,a,b);
    for(int i=0; i<l; ++i){
      t_n[i]=2.0*Ht_n_1[i]-t_n_2[i];
    }
    i_n*=i_minus;
    bessel=2.0*i_n*boost::math::cyl_bessel_j(n,a*t);
    for(int i=0; i<l; ++i){
      psi_t[i]+=bessel*t_n[i];
    }
    t_n_2=t_n_1;
    t_n_1=t_n;
  }
  
  faktor=-1.0i*(lambda_min+a*(1.0-epsilon/2.0))*t;
  faktor=exp(faktor);
  for(int i=0; i<l; ++i){
    psi_t[i]*=faktor;
  }
  
  return psi_t;
}
