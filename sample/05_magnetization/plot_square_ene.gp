set xlabel 'Magnetic Field'
set ylabel 'Energy'
set key left bottom
plot [0:5] 'energy_square.dat' u 1:2 w lp title "nstep = 100",\
 '' u 1:3 w lp title "nstep = 200",\
 '' u 1:4 w lp title "nstep = 500",\
 '' u 1:5 w lp title "nstep = 1000",\
 '' u 1:6 w lp title "nstep = 2000"
