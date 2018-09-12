#include "krylov.h"
#include "operatoren.h"

using namespace std;

vector<complex<double> > krylov(vector<complex<double> > psi,double t,int L, int m){
  int l=1<<L,k=0;
  vector<complex<double> > q_m,q_m1,phi;
  phi.resize(l);
  q_m.resize(l);
  q_m1.resize(l);
  complex<double> Q[l*m]={},d=0;
  double Norm=0.0,norm=0.0;
  //Variablen fuer zhetrd
  char uplo='U';
  complex<double> T[m*m]={};
  double D[m]={},E[m-1]={},TAU[m-1]={};
  complex<double> WORK_C[m*m]={};
  int info=0;
  //Variablen fuer dstev
  char v='V';
  double Z[m*m]={},WORK[2*m-2]={};
  
  //Normierung
  Norm=0;
  for(int i=0; i<l; ++i) Norm+=real(psi[i]*conj(psi[i]));
  Norm=sqrt(Norm);
  for(int i=0; i<l; ++i){
    Q[i*m]=psi[i]/Norm;
    psi[i]=psi[i]/Norm;
  }
  
  //Gram-Schmidt
  phi=H(psi,L);
  d=0;
  for(int i=0; i<l; ++i) d+=conj(psi[i])*phi[i];
  T[0]=d;
  norm=0;
  for(int i=0; i<l; ++i) norm+=real((phi[i]-d*psi[i])*conj(phi[i]-d*psi[i]));
  norm=sqrt(norm);
  for(int i=0; i<l; ++i) Q[i*m+1]=(phi[i]-d*psi[i])/norm;
  
  for(k=1; k<m-1; ++k){
    for(int i=0; i<l; ++i){
      q_m[i]=Q[i*m+k];
      q_m1[i]=Q[i*m+k-1];
    }
    
    phi=H(q_m,L);
    
    d=0;
    for(int i=0; i<l; ++i) d+=conj(q_m[i])*phi[i];
    T[k*m+k]=d;
    
    d=0;
    for(int i=0; i<l; ++i) d+=conj(q_m1[i])*phi[i];
    T[(k-1)*m+k]=d;
    T[k*m+k-1]=conj(d);
    norm=0;
    complex<double> normC=0;
    for(int i=0; i<l; ++i) normC+=(phi[i]-T[k*m+k]*q_m[i]-T[(k-1)*m+k]*q_m1[i])*
                          conj(phi[i]-T[k*m+k]*q_m[i]-T[(k-1)*m+k]*q_m1[i]);
    norm=real(normC);
    norm=sqrt(norm);
    for(int i=0; i<l; ++i) Q[i*m+k+1]=(phi[i]-T[k*m+k]*q_m[i]-T[(k-1)*m+k]*q_m1[i])/norm;
  }
  //T berechnen fuer k=m-1
  for(int i=0; i<l; ++i){
    q_m[i]=Q[i*m+k];
    q_m1[i]=Q[i*m+k-1];
  }
  phi=H(q_m,L);
  d=0;
  for(int i=0; i<l; ++i) d+=conj(q_m[i])*phi[i];
  T[k*m+k]=d;
  d=0;
  for(int i=0; i<l; ++i) d+=conj(q_m1[i])*phi[i];
  T[(k-1)*m+k]=d;
  T[k*m+k-1]=conj(d);
  
  //Komplexe Matrix in reelle Tridiagonalform umwandeln mit zhetrd
  zhetrd_(&uplo,&m,reinterpret_cast<double __complex__*>(&T[0]),&m,
	           reinterpret_cast<double*>(&D),
	           reinterpret_cast<double*>(&E),
	           reinterpret_cast<double __complex__*>(&TAU),
	           reinterpret_cast<double __complex__*>(&WORK_C),&m,&info);
  
  //Diese neue Matrix diagonalisieren mit dstev
  dstev_(&v,&m,D,E,Z,&m,WORK,&info);
  
  //Komplexe Matrizen
  complex<double> D_complex[m*m]={};
  complex<double> Z_complex[m*m]={};
  complex<double> Q_complex[l*m]={};
  complex<double> H[m*m]={};
  complex<double> C1[m*m]={},C2[m*m]={},C3[m*m]={},C4[l*m]={},alpha=1,beta=1;
  
  //Werte von D auf D_complex uebertragen und Exp berechnen
  for(int i=0; i<m;++i){
    D_complex[i*(m+1)]=-t*D[i]*1.0i;
    D_complex[i*(m+1)]=exp(D_complex[i*(m+1)]);
  }
  
  //Werte von Z auf Z_complex uebertragen
  for(int i=0; i<m*m;++i) Z_complex[i]=Z[i];
  
  //Werte von Q auf Q_complex uebertragen
  for(int i=0; i<l*m;++i) Q_complex[i]=Q[i];
  
  //Matrixmultiplikationen
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,m,m,
	      &alpha,&D_complex,m,&Z_complex,m,&beta,&C1,m);
  cblas_zgemm(CblasRowMajor,CblasConjTrans,CblasNoTrans,m,m,m,
	      &alpha,&Z_complex,m,&C1,m,&beta,&C2,m);
  
  for(int i=0; i<m*m; ++i) C1[i]=0;
  
  //Faktoren von TAU auf die Matrix anwenden
  for(int i=0; i<m-1; ++i) T[i*(m+1)+1]=1;
  
  for(k=m-2; k>-1; --k){
    for(int i=0; i<m; ++i) H[i*(m+1)]=1;
    for(int i=0; i<k+1; ++i){
      for(int j=0; j<k+1; ++j){
        H[i*m+j]-=TAU[k]*T[i*m+k+1]*conj(T[j*m+k+1]);
      }
    }
    cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,m,m,
	        &alpha,&C2,m,&H,m,&beta,&C1,m);
    cblas_zgemm(CblasRowMajor,CblasConjTrans,CblasNoTrans,m,m,m,
	        &alpha,&H,m,&C1,m,&beta,&C3,m);
    for(int i=0; i<m; ++i){
      H[i]=0;
      C2[i]=C3[i];
      C1[i]=0;
      C3[i]=0;
    }
  }
  
  //Multiplikation mit Q
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,l,m,m,
	      &alpha,&Q_complex,m,&C2,m,&beta,&C4,m);
  
  //Ausgabe
  vector<complex<double> > psi_aus;
  psi_aus.resize(l);
  for(int i=0; i<l; ++i) psi_aus[i]=Norm*C4[i*m];
  
  return psi_aus;
}
