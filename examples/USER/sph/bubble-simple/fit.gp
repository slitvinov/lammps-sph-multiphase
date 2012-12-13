set fit errorvariables
f(x) = ((x-x0)*(x>x0)*A)**B
A = 0.1
B = 1.5
C = 0.007
x0= 0.0318773177456613

set macros
t1=0.40
lim="[t1:][:]"
file20="data-nx20-ndim3-np1-D_heat_g0.05-alpha0.5-sph_c_g10.0-cg/rg.dat"
file40="data-nx40-ndim3-np1-D_heat_g0.05-alpha0.5-sph_c_g10.0-cg/rg.dat"
fit @lim f(x) file40 u 1:2 via A, x0

set xlabel "t, time"
set ylabel "Volume of bubble"


set key left
plot [:] file20 u 1:2 w l, \
     file40 u 1:2 w l, \
     f(x)/(x>t1) t sprintf("~ t^(%.2f)", B)
