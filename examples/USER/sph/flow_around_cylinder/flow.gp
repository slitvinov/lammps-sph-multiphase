set xlabel "y"
set ylabel "vx"
Lx=0.1
plot "<awk 'NF==4{print $2, $4}' data-wall/vcen.av | sort -g" u ($1*Lx):2, \
     "<awk 'NF==4{print $2, $4}' data-wall/vend.av | sort -g" u ($1*Lx):2, \
     "morris-fig6-1" w l, \
     "morris-fig6-2" w l 

