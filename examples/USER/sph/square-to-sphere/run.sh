#! /bin/bash

nx=20
ndim=3
np=2
dname=data-nx${nx}-ndim${ndim}-np${np}a
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp -var ndim ${ndim} -var nx ${nx} -var dname ${dname}
