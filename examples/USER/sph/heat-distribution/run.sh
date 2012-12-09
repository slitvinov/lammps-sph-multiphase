#! /bin/bash

nx=30
ndim=2
np=8
dname=data-nx${nx}-ndim${ndim}-np${np}
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp -var ndim ${ndim} -var nx ${nx} -var dname ${dname}

