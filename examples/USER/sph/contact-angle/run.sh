#! /bin/bash

nx=121
dname=data-wall-nx${nx}

# number of case
icase=1
mkdir -p ${dname}
mpirun -np 8  ../../../../src/lmp_linux -in droplet.lmp -var icase ${icase} -var dname ${dname} -var nx ${nx}