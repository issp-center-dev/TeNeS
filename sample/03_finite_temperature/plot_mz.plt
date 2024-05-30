set ytics 0.1
set xlabel "T" 
set ylabel "|<M_z>|" 

set xr [0.1:10.0]
set yr [-0.01:0.5]

p \
"magnetization_zero.dat" u (1/$1):(abs($2)) w l lc "black" t "h_x = 0",\
"magnetization_weak.dat" u (1/$1):(abs($2)) w l lc "red" t "h_x = 0.5",\
"magnetization_middle.dat" u (1/$1):(abs($2)) w l lc "green" t "h_x = 0.8",\
"magnetization_strong.dat" u (1/$1):(abs($2)) w l lc "blue" t "h_x = 2.0"
