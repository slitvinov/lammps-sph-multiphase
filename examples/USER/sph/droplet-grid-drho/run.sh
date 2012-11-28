#! /bin/bash

# number of processe
nproc=6
dname=data-wall
mkdir -p ${dname}
mpirun -np ${nproc}  ../../../../src/lmp_linux -in droplet.lmp -var dname ${dname}

