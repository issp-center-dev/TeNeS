set ytics 0.1
set xlabel "T" 
set ylabel "<M_x>" 

set log x

set xrange [0.1:10.0]
set yrange [-0.01:0.5]

p\
"magnetization_x_weak.dat" u (1/$1):2 w l lc "red" t "h_x = 0.5",\
"QMC/magx-L32-G0.5.dat" u (1/$1):2:3 w yer lc "red" pt 5 t "QMC h_x=0.5",\
"magnetization_x_middle.dat" u (1/$1):2 w l lc "green" t "h_x = 0.8",\
"QMC/magx-L32-G0.8.dat" u (1/$1):2:3 w yer lc "green"  pt 5 t "QMC h_x=0.8",\
"magnetization_x_strong.dat" u (1/$1):2 w l lc "blue" t "h_x = 2.0",\
"QMC/magx-L32-G2.0.dat" u (1/$1):2:3 w yer lc "blue" pt 5 t "QMC h_x=2.0"
