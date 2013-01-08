set style data l
system("./aux.sh")

dx(nx)=1.0/nx
R(V,nx) = (3.0/(4.0*pi))**(0.33333)*V**(0.3333)*dx(nx)
Rcor(V,nx) =  (3.0/(4.0*pi))**(0.33333)*V**(0.3333)*dx(nx) - dx(nx)

file="~/tmp/3d/data-nx40-ndim3-Lx1.0-D_heat_d0.6-alpha100-Hwv4.0-dprob0.01-time_k1.0-cv_d0.4-sph_rho_d0.01-dT0.5-cv_g0.04-D_heat_g0.1-sph_eta_d0.69395d/rg.dat"

set term x11 1
set xlabel "time"
set ylabel "<T_v>
plot "last1/rg.dat" u 1:($4/$6), \
     1.0, 0.0, 0.1

set term x11 2
set xlabel "time"
set ylabel "R^2"
plot [][:] \
     nx=40, "last1/rg.dat" u 1:(R($6,nx)**2), \
     "last1/rg.dat" u 1:(Rcor($6,nx)**2), \
     9.55*x - 1e-3, \
     1.4*(9.55*x - 1e-3)

set term x11 3
set xlabel "time"
set ylabel "R"
plot [][:] \
     nx=40, "last1/rg.dat" u 1:(R($6,nx)), \
     "last1/rg.dat" u 1:(Rcor($6,nx)), \
     sqrt(9.55*x)

     