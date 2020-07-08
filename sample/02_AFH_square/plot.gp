set key left top
set xlabel '1/D'
set ylabel '(E-Eqmc)/|Eqmc|'
set style arrow 1 size graph 0.02,20 linecolor rgb "black"
set dashtype 1 (5,5)
set arrow 1 nohead from 0,0.1666 to 5,0.1666 dt 1
exact = -0.6694422
plot [0:] \
 'energy.dat' u (1/$1):($2-exact)/-exact w lp title "FU = 0",\
 '' u (1/$1):($3-exact)/-exact w lp title "FU = 10",\
 '' u (1/$1):($4-exact)/-exact w lp title "FU = 100",\
# '' u (1/$1):($4-exact)/-exact w lp title "FU = 1000",\
