#! /bin/bash

nx=120
dname=data-wall-nx${nx}
mkdir -p ${dname}
mpirun -np 8  ../../../../src/lmp_linux -in droplet.lmp -var dname ${dname} -var nx ${nx}