#! /bin/bash
# x, y, fx, fy, colorgradient_peratom
../script/lammps2punto.sh data-wall/dump*.*  | awk 'NF{print $3, $4, $6, $7, $9, $10, $2, $NF} !NF' > hm
