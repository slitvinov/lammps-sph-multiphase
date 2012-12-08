#! /bin/bash

dname=data-wall3
mkdir -p ${dname}
mpirun -np 1  ../../../../src/lmp_linux -in insert.lmp -var dname ${dname}

