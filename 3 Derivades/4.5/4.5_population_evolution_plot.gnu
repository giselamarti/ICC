set title "Evolution of Population for Different k Values"
set xlabel "Iterations"
set ylabel "x"
set key top center outside maxrows 2
plot '4.5_population_evolution.dat' every ::2::22 u 1:2 w l lw 1.5 t "k = 0.5", \
     '4.5_population_evolution.dat' every ::25::45 u 1:2 w l lw 1.5 t "k = 1.0", \
     '4.5_population_evolution.dat' every ::48::68 u 1:2 w l lw 1.5 t "k = 1.5", \
     '4.5_population_evolution.dat' every ::71::91 u 1:2 w l lw 1.5 t "k = 2.0", \
     '4.5_population_evolution.dat' every ::94::114 u 1:2 w l lw 1.5 t "k = 2.5", \
     '4.5_population_evolution.dat' every ::117::137 u 1:2 w l lw 1.5 t "k = 3.0", \
     '4.5_population_evolution.dat' every ::140::160 u 1:2 w l lw 1.5 t "k = 3.5"