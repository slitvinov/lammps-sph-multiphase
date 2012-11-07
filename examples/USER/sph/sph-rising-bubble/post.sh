#! /bin/bash
# x, y, fx, fy, colorgradient_peratom
../script/lammps2punto.sh data/dump*.*  | awk 'NF{print $3, $4, $9, $10, $NF} !NF' > hm
