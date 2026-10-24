set size 0.8
set xlabel "time"
set ylabel "rel. shape anisotropy"

plot [:0.5][-0.0005:] \
     "data-wall-nx40/rg.dat" u 1:5 w l lw 4 t "nx=40",\
     "data-wall-nx68/rg.dat" u 1:5 w l lw 4 t "nx=69",\
     "data-wall-nx126/rg.dat" u 1:5 w l lw 4 t "nx=126"

call "../scripts/saver.gp" "cap_waves"