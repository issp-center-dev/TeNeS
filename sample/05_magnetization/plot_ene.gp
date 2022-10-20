set key left bottom
set xlabel 'Energy'
set ylabel 'Sz'
plot [0:5] 'energy.dat' u 1:2 w lp title "nstep = 100",\
 '' u 1:3 w lp title "nstep = 200",\
 '' u 1:4 w lp title "nstep = 500",\
 '' u 1:5 w lp title "nstep = 1000",\
 '' u 1:6 w lp title "nstep = 2000"
