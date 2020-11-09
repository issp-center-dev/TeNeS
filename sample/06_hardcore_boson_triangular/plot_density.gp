set key left top
set xlabel 'mu/V'
set ylabel 'n'
set y2label '|b|'
set style data lp

plot \
 'density.dat' u 1:4 pt 4 t"n",\
 'offdiag.dat' u 1:3 axes x1y2 pt 5 t"|b|"
