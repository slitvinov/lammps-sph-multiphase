#! /bin/bash

nx=40
np=8
ndim=3
D_heat_g=0.0005
alpha=0.50
sph_c_g=10.0

dname=data-nx${nx}-ndim${ndim}-np${np}-D_heat_g${D_heat_g}-alpha${alpha}-sph_c_g${sph_c_g}
rm -rf ${dname}
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp \
    -var alpha ${alpha} -var D_heat_g ${D_heat_g} -var ndim ${ndim} -var nx ${nx} \
    -var sph_c_g ${sph_c_g} -var dname ${dname}

