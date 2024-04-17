set multiplot layout 2,2
set title "Pointcare Map (k = 2.9)"
set xlabel "x(t)"
set ylabel "x(t+1)"
set xrange[0:1]
set yrange[0:1]
plot "4.5_poincare_map.dat" every ::2::22 u 1:2 t "k=2.9" w p pt 7 ps 0.7 lc rgb "black"
set title "Pointcare Map (k = 3.5)"
set xlabel "x(t)"
set ylabel "x(t+1)"
plot "4.5_poincare_map.dat" every ::25::45 u 1:2 t "k=3.5" w p pt 7 ps 0.7 lc rgb "black"
set title "Pointcare Map (k = 3.56)"
set xlabel "x(t)"
set ylabel "x(t+1)"
plot "4.5_poincare_map.dat" every ::48::68 u 1:2 t "k=3.56" w p pt 7 ps 0.7 lc rgb "black"
set title "Pointcare Map (k = 3.57)"
set xlabel "x(t)"
set ylabel "x(t+1)"
plot "4.5_poincare_map.dat" every ::71::91 u 1:2 t "k=3.57" w p pt 7 ps 0.7 lc rgb "black"
unset multiplot