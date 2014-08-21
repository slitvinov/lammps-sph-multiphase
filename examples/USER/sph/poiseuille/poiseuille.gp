
set xlabel "y"
set ylabel "vx"
set xtics 5e-4
plot "<awk 'NF==4{print $2, $4}' data/vx.av | sort -g" u ($1*2e-3):2 t "SPH", \
     "<seq 0 1e-5 1e-3 | ./poiseuille.awk -v t=5" w l t "theory"

call "../scripts/saver.gp" poiseuille