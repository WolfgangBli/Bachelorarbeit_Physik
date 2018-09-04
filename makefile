all: daten.out

daten.out: daten.cpp
	c++ -o daten.out -Wall daten.cpp -llapack -lblas

clean:
	rm daten.out
