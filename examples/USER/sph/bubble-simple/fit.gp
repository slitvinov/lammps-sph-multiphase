f(x) = A*abs(x)**B + C
A = 0.0038854976612842
B = 1.5
C = 1.0

set macros
lim="[0.06:0.15][:]"
file="data-nx20-ndim3-np8b/rg.dat"
fit @lim f(x) file u 1:2 via A, C, B

set xlabel "t, time"
set ylabel "Volume of bubble"

set fit errorvariables
set key left
plot file u 1:2 w lp, f(x) t sprintf("A*t^%.2f +/- %.3f", B, B_err)