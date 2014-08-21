set macro
dt = system("cut -d ' ' -f 4 in.dt") 
set term x11 1
plot [:0.5][-0.5:0.5] "vom.dat" u (dt*$1):2 w lp, "" u (dt*$1):3 w lp

set term x11 2
case1="[0.065:0.095]"
case2=""
plot [0:0.5]@case2 "com.dat" u (dt*$1):2 w lp, "" u (dt*$1):3 w lp

