#! /bin/bash

gy=$1
dname=data-wall-gy${gy}-part3
mkdir -p ${dname}
mpirun -np 1  ../../../../src/lmp_linux -in insert.lmp -var gy ${gy} -var dname ${dname}

