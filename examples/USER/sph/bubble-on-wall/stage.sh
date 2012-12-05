#! /bin/bash

for g in 1.0 1.25 1.5 2.0 3.0; do
    Np=$(awk '$3>1.0&&NR>1{print $2;exit}' data-wall-gy${g}/rg.dat)
    echo $g ${Np}
done


