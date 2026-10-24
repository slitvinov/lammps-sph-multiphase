set macro
dt = system("cut -d ' ' -f 4 in.dt") 
set term x11 1
plot "data-ndim2/rg.dat" w lp
