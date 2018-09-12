all: daten.out

daten.out: daten.o operatoren.o krylov.o chebyshev.o
	c++ -Wall -o daten.out daten.o operatoren.o krylov.o chebyshev.o -llapack -lblas

daten.o: daten.cpp
	g++ -Wall -c daten.cpp -llapack -lblas

operatoren.o: operatoren.cpp
	g++ -Wall -c operatoren.cpp

krylov.o: krylov.cpp
	g++ -Wall -c krylov.cpp -llapack -lblas

chebyshev.o: chebyshev.cpp
	g++ -Wall -c chebyshev.cpp

clean:
	rm daten.out
	rm daten.o
	rm operatoren.o
	rm krylov.o
	rm chebyshev.o
