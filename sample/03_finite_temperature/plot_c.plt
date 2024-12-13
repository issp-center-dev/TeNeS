set auto
unset log
set log x

set xr [0.1:]
set yr [0:]

set xl "temperature"
set yl "specific heat"

plot \
"spec_zero.dat" u (1/$1):2 w l lc "black" t"iTPS h_x=0.0", \
"QMC/spec-L32-G0.0.dat" u (1/$1):2 w yer lc "black" pt 5 t"QMC h_x=0.0", \
"spec_weak.dat" u (1/$1):2 w l lc "red" t"iTPS h_x=0.5", \
"QMC/spec-L32-G0.5.dat" u (1/$1):2 w yer lc "red" pt 5 t"QMC h_x=0.5", \
"spec_middle.dat" u (1/$1):2 w l lc "green" t"iTPS h_x=0.8", \
"QMC/spec-L32-G0.8.dat" u (1/$1):2 w yer lc "green" pt 5 t"QMC h_x=0.8", \
"spec_strong.dat" u (1/$1):2 w l lc "blue" t"iTPS h_x=2.0", \
"QMC/spec-L32-G2.0.dat" u (1/$1):2 w yer lc "blue" pt 5 t"QMC h_x=2.0", \
