#! /bin/bash

nx=40
np=8
ndim=3
D_heat_g=0.05
cp=0.50
alpha=0.5
sph_c_g=10.0
prob=1.0
time_k=4.00

dname=data-nx${nx}-ndim${ndim}-np${np}-D_heat_g${D_heat_g}-alpha${alpha}-sph_c_g${sph_c_g}-cp${cp}-prob${prob}-time_k${time_k}
rm -rf ${dname}
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp \
    -var alpha ${alpha} -var D_heat_g ${D_heat_g} -var ndim ${ndim} -var nx ${nx} \
    -var sph_c_g ${sph_c_g} -var cp ${cp} -var prob ${prob} -var time_k ${time_k} \
    -var dname ${dname}

