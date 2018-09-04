set terminal pdf color size 21cm,29.7cm
set output "Zeitentwicklung_S_z(0).pdf"
set style data linespoints

plot "data/krylov_L11_m5_t0.1_T20.data" using 1:3 title "S_z(0) in Abängigkeit von t"
