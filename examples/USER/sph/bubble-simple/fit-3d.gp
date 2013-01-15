set macro
dx(nx)=1.0/nx
R(V,nx) = (3.0/(4.0*pi))**(0.33333)*V**(0.3333)*dx(nx)
st = '1:( R($6, nx)**2 )

set style data l
set key left
set term x11 1

set xlabel "t"
set ylabel "R^2"

set size 0.75
plot [:2][0:] \
     nx=20, "/scratch/data/3d-bubble-wall/data-nx20-ndim3-Lx1.0-D_heat_g0.04-alpha0.25-Hwv1.0-dprob0.01-time_k1.0-cv_g4.0-sph_rho_g100-dT0.01/rg.dat" u @st w l lw 3 t "nx=20", \
     nx=40, "/scratch/data/3d-bubble-wall/data-nx40-ndim3-Lx1.0-D_heat_g0.04-alpha0.25-Hwv1.0-dprob0.01-time_k1.0-cv_g4.0-sph_rho_g100-dT0.01/rg.dat" u @st w l lw 3 t "nx=40", \
     nx=80, "/scratch/data/3d-bubble-wall/data-nx80-ndim3-Lx1.0-D_heat_g0.04-alpha0.25-Hwv1.0-dprob0.01-time_k1.0-cv_g4.0-sph_rho_g100-dT0.01/rg.dat" u @st w l lw 3 t "nx=80", \
     4e-2*x - 0.01 t "Mikic et al. (1970)"

     
call "saver.gp" "f1"

# set term x11 2
# plot "/scratch/data/3d-bubble-wall/data-nx20-ndim3-Lx1.0-D_heat_g0.04-alpha0.25-Hwv1.0-dprob0.01-time_k1.0-cv_g4.0-sph_rho_g100-dT0.01/rg.dat" w l lw 3 t "nx=20", \
#      "/scratch/data/3d-bubble-wall/data-nx40-ndim3-Lx1.0-D_heat_g0.04-alpha0.25-Hwv1.0-dprob0.01-time_k1.0-cv_g4.0-sph_rho_g100-dT0.01/rg.dat" w l lw 3 t "nx=40", \
#      "/scratch/data/3d-bubble-wall/data-nx80-ndim3-Lx1.0-D_heat_g0.04-alpha0.25-Hwv1.0-dprob0.01-time_k1.0-cv_g4.0-sph_rho_g100-dT0.01/rg.dat" w l lw 3 t "nx=80"

# call "saver.gp" f2