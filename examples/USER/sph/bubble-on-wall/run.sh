#! /bin/bash

gy=$1
dname=data-wall-gy${gy}
mkdir -p ${dname}
mpirun -np 8  ../../../../src/lmp_linux -in insert.lmp -var gy ${gy} -var dname ${dname}

