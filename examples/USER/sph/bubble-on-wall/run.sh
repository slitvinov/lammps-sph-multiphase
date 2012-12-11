#! /bin/bash

nx=$1
gy=$2
np=8
dname=data-wall-gy${gy}-nx${nx}-np${np}
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp -var nx ${nx} -var gy ${gy} -var dname ${dname}

