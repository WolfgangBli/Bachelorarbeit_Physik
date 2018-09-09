# Bachelorarbeit Physik
Dies ist eine Zusammenfassung meiner Bachelorarbeit in Physik in der ich ein Programm zur Simulation einer Spinkette geschrieben habe.  
Die Datei [daten.cpp](https://github.com/WolfgangBli/Bachelorarbeit_Physik/blob/master/daten.cpp) enthält den Quelltext. Für diese Version wurden kleinere Vereinfachungen vorgenommen um die wesentlichen Funktionen hervorzuheben. Zur besseren Verständlichkeit der Definitionen der Funktionen *H* und *S_z* im Quelltext sind die [physikalischen Grundlagen](https://github.com/WolfgangBli/Bachelorarbeit_Physik/blob/master/Physikalische_Grundlagen.pdf "Physikalische_Grundlagen.pdf") kurz zuasammengefasst. Wer sich auch für die numerischen Verfahren auf denen die Funktionen *krylov* und *chebyshev* beruhen interessiert findet diese in meiner [vollständigen Arbeit](https://github.com/WolfgangBli/Bachelorarbeit_Physik/blob/master/Vollstaendige_Arbeit.pdf "Vollstaendige_Arbeit.pdf"). Das Programm erzeugt Daten, die im Verzeichnis *data* gespeichert werden. Für die Arbeit wurden aus diesen Daten Graphiken erzeugt. Als Beispiel wurde hier eine solche [graphische Auswertung](https://github.com/WolfgangBli/Bachelorarbeit_Physik/blob/master/Zeitentwicklung_S_z0.pdf "Zeitentwicklung_S_z0.pdf") hinzugefügt, welche mit *gnuplot* und [Zeitentwicklung.gp](https://github.com/WolfgangBli/Bachelorarbeit_Physik/blob/master/Zeitentwicklung.gp) erzeugt wurde.

## Programmaufruf
Beim Programmaufruf werden Parameter des Systems und der numerischen Verfahren übergeben. Ein Programmaufruf sieht wie folgt aus:

    ./daten.out L m N t_k t_c T Methode #

Die Bedeutung der einzelnen Variablen ist in der Liste unten aufgeführt. Ein typischer Programmaufruf sieht zum Beispiel so aus:

    ./daten.out 11 5 5 0.1 0.1 20 k 1

Dieser Aufruf führt zu einer Programmlaufzeit von einigen Sekunden. **Vorsicht**, bei höheren Werten von *L*, *m* oder *N* kann die Programmlaufzeit sehr viel länger sein.

| Variable | Datentyp | Bedeutung                                             |
| :------- | :------- | :---------------------------------------------------- |
| L        | int      | Anzahl der Spins auf der Spinkette                    |
| m        | int      | Parameter Krylov                                      |
| N        | int      | Parameter Chebyshev                                   |
| t_k      | double   | Schrittweite Krylov                                   |
| t_c      | double   | Schrittweite Chebyshev                                |
| T        | double   | Gesamtzeit die zu simmulieren ist                     |
| Methode  | char     | Verwendete Methode, mögliche Werte: 'k', 'c' oder 'b' |
| #        | int      | Zahl zur Nummerierung in der Ausgabedatei             |
