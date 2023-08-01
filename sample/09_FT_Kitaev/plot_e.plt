set auto
unset log
set log x
set key left top

set xl 'temperature'
set yl 'energy density'

set xr [0.01:20]

plot \
'TPQ/TPQ_gapless.dat' u (1/$1):2 w l lw 2 lc rgb 'black' t 'gapless, TPQ', \
'energy_gapless.dat' u (1/$1):2 pt 5 lc rgb 'blue' t 'gapless, iTPS', \
'TPQ/TPQ_gapfull.dat' u (1/$1):2 w l lw 2 lc rgb 'black' dt (10,5) t 'gapfull, TPQ', \
'energy_gapfull.dat' u (1/$1):2 pt 7 lc rgb 'red' t 'gapfull, iTPS', \
