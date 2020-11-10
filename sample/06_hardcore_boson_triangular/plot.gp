set key left top
set xlabel 'mu/V'
set ylabel 'S(Q)/N'
set y2label '|b|'
set style data lp

plot \
 'sq.dat' u 1:3 pt 4 t"S(Q)",\
 'offdiag.dat' u 1:3 axes x1y2 pt 5 t"|b|"
