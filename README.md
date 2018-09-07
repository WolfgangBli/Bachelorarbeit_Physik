# Bachelorarbeit Physik
Dies ist der Quelltext für das Programm meiner Bachelorarbeit in Physik.
Ein Programmaufruf sieht wie folgt aus:

    ./daten.out L m N t_k t_c T Methode #

Dabei sind

| Variable | Datentyp | Bedeutung              |
| :------- | :------- | :--------------------- |
| L        | int      | Anzahl Spins           |
| m        | int      | Parameter Krylov       |
| N        | int      | Parameter Chebyshev    |
| t_k      | double   | Schrittweite Krylov    |
| t_c      | double   | Schrittweite Chebyshev |
| T        | double   | Gesamtzeit             |
| Methode  | char     | Verwendete Methode     |
| #        | int      | Zahl zum nummerieren   |
