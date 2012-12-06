#! /bin/bash

../script/lammps2punto.sh $1/dump* | awk 'NF{print $3, $4, $5, $6, $7, $8, $2, $NF} !NF' > $1/temp.dat
