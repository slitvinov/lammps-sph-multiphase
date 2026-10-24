set xlabel "y"
set ylabel "vx"
Lx=0.1
set key center
plot [0:Lx] \
     "<awk 'NF==4{print $2, $4}' data-wall/vcen.av | sort -g" u ($1*Lx):2 t "SPH", \
     "<awk 'NF==4{print $2, $4}' data-wall/vend.av | sort -g" u ($1*Lx):2 t "", \
     "morris-fig6-1" w l t "FEM, Path 1", \
     "morris-fig6-2" w l t "FEM, Path 2"

call "../scripts/saver.gp" flow
