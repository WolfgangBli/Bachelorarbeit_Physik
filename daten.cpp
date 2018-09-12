#include "chebyshev.h"
#include "operatoren.h"
#include "krylov.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>

#include <time.h>

using namespace std;

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
