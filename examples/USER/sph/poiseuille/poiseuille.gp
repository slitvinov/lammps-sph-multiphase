set size 0.8

set xlabel "poistion"
set ylabel "vx, velocity"
set xtics 1e-3
plot "<awk 'NF==4{print $2, $4}' data/vx.av | sort -g" u ($1*2e-3):2 w p ps 3 t "SPH", \
     "<seq 0 1e-5 1e-3 | ./poiseuille.awk -v t=5" w l lw 4 t "theory"

call "../scripts/saver.gp" poiseuille