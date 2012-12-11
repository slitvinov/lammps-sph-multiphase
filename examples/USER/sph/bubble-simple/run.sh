#! /bin/bash

nx=$1
ndim=3
np=8
dname=data-nx${nx}-ndim${ndim}-np${np}-wall
rm -rf ${dname}
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp -var ndim ${ndim} -var nx ${nx} -var dname ${dname}

