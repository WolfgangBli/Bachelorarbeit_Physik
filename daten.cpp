#include <boost/math/special_functions/bessel.hpp> 
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
#include <lapacke.h>
#include <cblas.h>
#include <fstream>
#include <sstream>
#include <time.h>

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
	default: cout<<"Fehler"<<endl; break;
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
      default: cout<<"Fehler"<<endl; break;
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


int main(int argn, char* argv[]){
  /*Physikalische Kenngroessen des Systems:
   * L ist die Anzahl der Elektronenspins auf der Kette.
   * l ist die Anzahl Basisvektoren im zugehoerigen Hilbertraum,
     es gilt l=2^L.*/
  int L,l;
  /*Parameter für die numerischen Verfahren,
    je groesser diese Werte, desto genauer sind die Verfahren:
   * m ist die Dimension des Krylovraums.
   * N legt die fest, bis zu welchem Punkt die Partialsumme beim 
     Chebyshevverfahren berechnet wird.*/
  int m,N;
  /*number dient lediglich dazu, um die einzelnen Programmaufrufe in der
    Ausgabedatei Zeiten.data nummerieren zu koennen.*/
  int number;
  /*Zeitparameter des Systems. Ziel der Programms ist es die Zeitentwicklung
    des Systems zu simmulieren, dies wird durch iteratives anwenden der
    Verfahren erreicht.
   * T ist die absolute Zeit, bis zu der die Zeitentwicklung bestimmt wird,
     typischer Wert ist 20.
   * t_k, bzw. t_c ist die Groesse der Iterationsschritte fuer das Krylov- bzw.
     Chebyshevverfahren, typischer Wert ist 0.1.*/
  double t_k,t_c,T;
  /*Mit meth wird die numerische Methode festgelegt,
    die Verwendet werden soll.
    'k' fuer Krylov,
    'c' fuer Chebyshev oder
    'b' fuer beide.*/
  char meth;
  /*Anzahl der benoetigten Iterationsschritte fuer die jeweiligen Verfahren,
    berechnet sich aus T und t_k bzw. t_c*/
  int Zeitschritte_k,Zeitschritte_c;
  //Parameter einlesen
  if (argn!=9){
    cout<<"Eingabe: L m N t_k t_c T Methode #"<<endl;
    return 0;
  }
  istringstream ss1(argv[1]);
  if (!(ss1>>L)) cerr << "Invalid number " << argv[1] << '\n';
  istringstream ss2(argv[2]);
  if (!(ss2>>m)) cerr << "Invalid number " << argv[1] << '\n';
  istringstream ss3(argv[3]);
  if (!(ss3>>N)) cerr << "Invalid number " << argv[1] << '\n';
  istringstream ss4(argv[4]);
  if (!(ss4>>t_k)) cerr << "Invalid number " << argv[1] << '\n';
  istringstream ss5(argv[5]);
  if (!(ss5>>t_c)) cerr << "Invalid number " << argv[1] << '\n';
  istringstream ss6(argv[6]);
  if (!(ss6>>T)) cerr << "Invalid number " << argv[1] << '\n';
  istringstream ss7(argv[7]);
  if (!(ss7>>meth)) cerr << "Invalid char " << argv[1] << '\n';
  istringstream ss8(argv[8]);
  if (!(ss8>>number)) cerr << "Invalid number " << argv[1] << '\n';
  
  l=1<<L;
  Zeitschritte_k=T/t_k+1;
  Zeitschritte_c=T/t_c+1;
  
  /*betrag ist eine positve reele Zahl, die zur Messung verschiedener
    Observablen verwendet wird. Sie wird als complex<double> gespeichert um mit
    den folgenden Vektoren kompatibel zu sein. Zur Ausgabe wird nur der Realteil
    verwendet.*/
  complex<double> betrag=0;
  
  /*In den Vektoren psi_k und psi_c wird der Zustand des physikalischen Systems
    gespeichert. phi dient als Zwischenspeicher.*/
  vector<complex<double> > psi_k,psi_c,phi;
  phi.resize(l);
  psi_k.resize(l);
  psi_c.resize(l);
  /*Als Initialisierung wird hier ein Neel-Zustand gesetzt.*/
  int index=0;
  for(int k=1; k<L; k+=2){
    index+=1<<k;
  }
  psi_k[index]=1;
  psi_c[index]=1;
  
  /*Variablen fuer die Messung der Programmlaufzeit.*/
  double time_k=0.0,time_c=0.0,tstart;
  
  
  //Textausgabe
  stringstream dateiname_k,dateiname_c;
  dateiname_k<<"data/krylov_L"<<L<<"_m"<<m<<"_t"<<t_k<<"_T"<<T<<".data";
  dateiname_c<<"data/chebyshev_L"<<L<<"_N"<<N<<"_t"<<t_c<<"_T"<<T<<".data";
  ofstream fileout;
  
  //Krylov aufrufen
  if(meth=='k'||meth=='b'){
    tstart=clock();
    fileout.open((dateiname_k.str()).c_str());
    fileout<<"#Krylov L = "<<L<<" m = "<<m<<endl;
    for(int k=0; k<Zeitschritte_k; k++){
      //Energie <psi|H|psi>
      phi=H(psi_k,L);
      betrag=0;
      for(int i=0; i<l; i++){
        betrag+=conj(psi_k[i])*phi[i];
      }
      fileout<<t_k*k<<"   "<<real(betrag);
      //<psi|S_z(j)|psi>
      for(int j=0; j<L; j++){
        phi=S_z(psi_k,L,j);
        betrag=0;
        for(int i=0; i<l; i++){
          betrag+=conj(psi_k[i])*phi[i];
        }
        fileout<<"   "<<real(betrag);
      }
      fileout<<endl;
      psi_k=krylov(psi_k,t_k,L,m);
    }
    fileout.close();
    time_k=clock()-tstart;
    time_k=time_k/CLOCKS_PER_SEC;
    cout<<"Zeit Krylov:"<<time_k<<endl;
  }
  
  //Chebyshev aufrufen
  if(meth=='c'||meth=='b'){
    tstart=clock();
    fileout.open((dateiname_c.str()).c_str());
    fileout<<"#Chebyshev L = "<<L<<" N = "<<N<<endl;
    for(int k=0; k<Zeitschritte_c; k++){
      //Energie <psi|H|psi>
      phi=H(psi_c,L);
      betrag=0;
      for(int i=0; i<l; i++){
        betrag+=conj(psi_c[i])*phi[i];
      }
      fileout<<t_c*k<<"   "<<real(betrag);
      //<psi|S_z(j)|psi>
      for(int j=0; j<L; j++){
        phi=S_z(psi_c,L,j);
        betrag=0;
        for(int i=0; i<l; i++){
          betrag+=conj(psi_c[i])*phi[i];
        }
        fileout<<"   "<<real(betrag);
      }
      fileout<<endl;
      psi_c=chebyshev(psi_c,t_c,L,N);
    }
    fileout.close();
    time_c=clock()-tstart;
    time_c=time_c/CLOCKS_PER_SEC;
    cout<<"Zeit Chebyshev:"<<time_c<<endl;
  }
  
  return 0;
}
