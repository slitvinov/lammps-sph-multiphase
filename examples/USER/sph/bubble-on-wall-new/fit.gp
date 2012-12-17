set fit errorvariables
f(x) = ((x-x0)*(x>x0)*A)**B + C
A = 1.05867741953623
B = 1.5
C = 0.007
x0= 0.0642258364284927

set macros
t1=0.2
lim="[t1:]"
file20="data-nx20-ndim3-np1-D_heat_g0.05-alpha0.5-sph_c_g10.0newplace/rg.dat"
file40="data-nx40-ndim3-np1-D_heat_g0.05-alpha0.5-sph_c_g10.0newplace/rg.dat"
fit @lim f(x) file40 via C, B, x0, A

set xlabel "t, time"
set ylabel "Volume of bubble"


set key left
plot [:] file20 u 1:2 w l lw 3, \
     file40 u 1:2 w l lw 3, \
     f(x) t sprintf("~ t^(%.2f +/ %.3f)", B, B_err)
