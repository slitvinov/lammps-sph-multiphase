#! /bin/bash

dname=data-wall
mkdir -p ${dname}
mpirun -np 8  ../../../../src/lmp_linux -in poiseuille.lmp -var dname ${dname}

