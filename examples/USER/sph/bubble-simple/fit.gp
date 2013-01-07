system("./aux.sh")

dx(nx)=1.0/nx
R(V,nx) = (3.0/(4.0*pi))**(0.33333)*V**(0.3333)*dx(nx)
Rcor(V,nx) =  (3.0/(4.0*pi))**(0.33333)*V**(0.3333)*dx(nx) - dx(nx)

file40="~/tmp/3d/data-nx40-ndim3-Lx1.0-D_heat_d0.6-alpha100-Hwv4.0-dprob0.01-time_k1.0-cv_d1.0-sph_rho_d0.01-dT0.1-cv_g0.04-D_heat_g0.1-sph_eta_d0.69395/rg.dat"
file80="~/tmp/3d/data-nx80-ndim3-Lx1.0-D_heat_d0.6-alpha100-Hwv4.0-dprob0.01-time_k1.0-cv_d1.0-sph_rho_d0.01-dT0.1-cv_g0.04-D_heat_g0.1-sph_eta_d0.69395/rg.dat"

set style data l
plot [][:] \
     nx=20, "last1/rg.dat" u 1:(R($6,nx)**2), \
     (9.549296585513725*x-0.01)*1.0

     #plot "last1/rg.dat" u ($3/$6), 1.0, 1.5