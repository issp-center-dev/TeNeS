set ytics 0.1
set xlabel "T" 
set ylabel "<M_x>" 

p [:5] \
"magnetization_x_zero.dat" u (1/$1):2 w lp ti "h_x = 0",\
"magnetization_x_weak.dat" u (1/$1):2 w lp ti "h_x = 0.5",\
"magnetization_x_middle.dat" u (1/$1):2 w lp ti "h_x = 0.8",\
"magnetization_x_strong.dat" u (1/$1):2 w lp ti "h_x = 2.0"
