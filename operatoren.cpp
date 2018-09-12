#include "operatoren.h"

using namespace std;

vector<complex<double> > H(vector<complex<double> > psi,int L){
  int j=0,l=1<<L;
  vector<complex<double> > phi;
  phi.resize(l);
  //Schleife fuer alle Summanden des Hamilton
  for(int k=0; k<L-1; ++k){
    //Schleife fuer die Vektorkomponenten
    for(int i=0; i<l; ++i){
      switch((i>>k)&3){
          case 0: phi[i]+=0.25*psi[i];
            break;
          case 1: j = i+(1<<k);
            phi[j]+=0.5*psi[i];
            phi[i]+=(-0.25)*psi[i];
            break;
          case 2: j = i-(1<<k);
            phi[j]+=0.5*psi[i];
            phi[i]+=(-0.25)*psi[i];
            break;
          case 3: phi[i]+=0.25*psi[i];
            break;
          default: break;
      }
    }
  }
  //Periodische Randbedingungen fuer gerades L
  if((L&1)==0){
    //Schleife fuer die Vektorkomponenten
    for(int i=0; i<l; ++i){
      if(((i&1)==0)&&(((i>>(L-1))&1)==0)) phi[i]+=0.25*psi[i];
      if(((i&1)==0)&&(((i>>(L-1))&1)==1)){
        j=i-(l>>1)+1;
        phi[j]+=0.5*psi[i];
        phi[i]+=(-0.25)*psi[i];
      }
      if(((i&1)==1)&&(((i>>(L-1))&1)==0)){
        j=i+(l>>1)-1;
        phi[j]+=0.5*psi[i];
        phi[i]+=(-0.25)*psi[i];
      }
      if(((i&1)==1)&&(((i>>(L-1))&1)==1)) phi[i]+=0.25*psi[i];
    }
  }
  return phi;
}

vector<complex<double> > S_z(vector<complex<double> > psi,int L,int i){
  int l=1<<L;
  vector<complex<double> > phi;
  phi.resize(l);
  //Schleife fuer die Vektorkomponenten
  for(int k=0; k<l; ++k){
    switch((k>>i)&1){
      case 0: phi[k]+=(-0.5)*psi[k];
        break;
      case 1: phi[k]+=0.5*psi[k];
        break;
      default: break;
    }
  }
  return phi;
}

vector<complex<double> > H_rescaled(vector<complex<double> > psi,int L,double a, double b){
  int l=pow(2,L);
  vector<complex<double> > phi;
  phi.resize(l);
  phi=H(psi,L);
  for(int i=0; i<l; ++i){
    phi[i]=(phi[i]-b*psi[i])/a;
  }
  return phi;
}
