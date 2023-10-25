set auto
unset log
set log x

set xr [0.1:]

set xl "temperature"
set yl "energy density"

plot \
"energy_zero.dat" u (1/$1):2 w l lc "black" t"iTPS h_x=0.0", \
"QMC/ene-L32-Gamma0.0.dat" w yer lc "black" pt 5 t"QMC h_x=0.0", \
"energy_weak.dat" u (1/$1):2 w l lc "red" t"iTPS h_x=0.5", \
"QMC/ene-L32-Gamma0.5.dat" w yer lc "red" pt 5 t"QMC h_x=0.5", \
"energy_middle.dat" u (1/$1):2 w l lc "green" t"iTPS h_x=0.8", \
"QMC/ene-L32-Gamma0.8.dat" w yer lc "green" pt 5 t"QMC h_x=0.8", \
"energy_strong.dat" u (1/$1):2 w l lc "blue" t"iTPS h_x=2.0", \
"QMC/ene-L32-Gamma2.0.dat" w yer lc "blue" pt 5 t"QMC h_x=2.0", \
