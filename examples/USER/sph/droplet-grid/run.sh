#! /bin/bash

dname=data-wall
mkdir -p ${dname}
mpirun -np 6  ../../../../src/lmp_linux -in droplet.lmp -var dname ${dname}

