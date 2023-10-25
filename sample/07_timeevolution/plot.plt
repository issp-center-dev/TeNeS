set yr [-0.5:0.5]
set ytics 0.25
set key bottom
set xl "time"
set yl "magnetization"

plot \
"magnetization_strong.dat" w lp pt 5 t"hx=2.0", \
"magnetization_middle.dat" w lp pt 7 t"hx=0.8", \
"magnetization_weak.dat" w lp pt 9 t"hx=0.5", \
0.0 lc "black" t""
