set multiplot layout 3,3
unset key
set xlabel "t"

set title "Gaussian function (signal)"
set ylabel "f(t)"
plot "funcions.dat" u 1:2 w lp pt 7 ps 0.5 lc rgb "black"

set title "Noise"
set ylabel "noise(t)"
plot "funcions.dat" u 1:3 w lp pt 7 ps 0.5 lc rgb "black"

set title "Signal and noise"
set ylabel "singal+noise(t)"
plot "funcions.dat" u 1:4 w lp pt 7 ps 0.5 lc rgb "black"

set title "FT of the Gaussian function (signal)"
set ylabel "f(t)"
plot "FT.dat" u 1:2 w lp pt 7 ps 0.5 lc rgb "black"

set title "FT of the Noise"
set ylabel "noise(t)"
plot "FT.dat" u 1:3 w lp pt 7 ps 0.5 lc rgb "black"

set title "FT of the Signal and noise"
set ylabel "singal+noise(t)"
plot "FT.dat" u 1:4 w lp pt 7 ps 0.5 lc rgb "black"

set title "IFT of the Gaussian function (signal)"
set ylabel "f(t)"
plot "InvFT.dat" u 1:2 w lp pt 7 ps 0.5 lc rgb "black"

set title "IFT of the Noise"
set ylabel "noise(t)"
plot "InvFT.dat" u 1:3 w lp pt 7 ps 0.5 lc rgb "black"

set title "IFT of the Signal and noise"
set ylabel "singal+noise(t)"
plot "InvFT.dat" u 1:4 w lp pt 7 ps 0.5 lc rgb "black"