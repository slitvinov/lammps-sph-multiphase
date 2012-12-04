#! /bin/bash

../script/lammps2punto.sh data-wall/dump* | awk 'NF{print $3, $4, $5, $NF} !NF' > temp.dat
