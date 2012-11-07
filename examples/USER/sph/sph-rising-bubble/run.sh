#! /bin/bash

dname=data-wall
mkdir -p ${dname}
mpirun -np 2  ../../../../src/lmp_linux -in bubble.lmp -var dname ${dname}

