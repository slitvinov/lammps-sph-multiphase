set macro
dt = 1.7146336010054037343e-06
set term x11 1
plot [:0.5][-0.5:0.5] "vom.dat" u (dt*$1):2 w lp, "" u (dt*$1):3 w lp

set term x11 2
case1="[0.065:0.095]"
#case2="[0.081:0.087]"
case2=""
plot [0:0.5]@case2 "com.dat" u (dt*$1):2 w lp, "" u (dt*$1):3 w lp

