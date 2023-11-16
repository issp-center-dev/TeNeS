set key left top
# set xlabel "mu/V"
# set ylabel "S(Q), N(Q)"
# set y2label "|b|+|b'|"
set ytics nomirror
set y2tics nomirror
set style data lp

plot \
 'nq.dat' u 1:(abs($3)) pt 4 t"|N(Q)|",\
 'sq.dat' u 1:(abs($3)) pt 5 t"|S(Q)|",\
 'sq_x.dat' u 1:(abs($3)) pt 6 t"|S'(Q)|",\
 'offdiag.dat' u 1:(abs($3)) axes x1y2 pt 7 t"|b|+|b'|"
