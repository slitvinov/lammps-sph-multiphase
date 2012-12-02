#! /bin/bash

dname=data-wall
mkdir -p ${dname}
mpirun -np 1  ../../../../src/lmp_linux -in poiseuille.lmp -var dname ${dname}

