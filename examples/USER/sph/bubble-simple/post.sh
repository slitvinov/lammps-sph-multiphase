#! /bin/bash

../script/lammps2punto.sh $1/dump* | awk 'NF{print $3, $4, $5, $6, $7, $9, $NF, $2} !NF' > $1/temp.dat
../script/lammps2punto.sh $1/dump* | awk 'NF&&$2==2{print $3, $4, $5, $6, $7, $9, $NF, $2} !NF' > $1/bubble.dat
../script/lammps2punto.sh $1/dump* | awk 'NF&&$3<0.5{print $3, $4, $5, $6, $7, $9, $NF, $2} !NF' > $1/half.dat

../script/lammps2punto.sh $1/dump* | awk 'NF&&$3<0.80&&$3>0.70{print $3, $4, $5, $6, $7, $9, $NF, $2} !NF' > $1/slice.x.dat
../script/lammps2punto.sh $1/dump* | awk 'NF&&$4<0.55&&$4>0.45{print $3, $4, $5, $6, $7, $9, $NF, $2} !NF' > $1/slice.y.dat
../script/lammps2punto.sh $1/dump* | awk 'NF&&$5<0.55&&$5>0.45{print $3, $4, $5, $6, $7, $9, $NF, $2} !NF' > $1/slice.z.dat

awk 'f{print $1; f=0; nextfile} /ITEM: NUMBER OF ATOMS/{f=1}' $1/dump* > $1/n.dat
