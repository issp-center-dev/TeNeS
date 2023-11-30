set ytics 0.1
set xlabel "T" 
set ylabel "|<M_z>|" 

p [:5] \
"magnetization_zero.dat" u (1/$1):(abs($2)) w lp ti "h_x = 0",\
"magnetization_weak.dat" u (1/$1):(abs($2)) w lp ti "h_x = 0.5",\
"magnetization_middle.dat" u (1/$1):(abs($2)) w lp ti "h_x = 0.8",\
"magnetization_strong.dat" u (1/$1):(abs($2)) w lp ti "h_x = 2.0"
