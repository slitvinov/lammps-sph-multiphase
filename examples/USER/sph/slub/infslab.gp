set xlabel "position"
set ylabel "temerature"

plot "infslab.dat" w l t "reference", \
     "dump.last" u 3:6 w p t "SPH"

call "../scripts/saver.gp" infslab     