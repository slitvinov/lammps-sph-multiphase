#! /bin/bash

Rc=0.002
function getnp_depart() {
    awk '$8>1.0&&NR>1{print $2;exit}' $1
}

function getg() {
    echo  $1 | sed -e 's/.*gy//1' -e 's/-.*//1'
}


for d in $(ls -d /scratch/data/2dlong/data*); do
    g=$(getg $d)
    Np=$(getnp_depart $d/rg.dat)
    echo $g ${Np}
done | awk 'NF==2'


