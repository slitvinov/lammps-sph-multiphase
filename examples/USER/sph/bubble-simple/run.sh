#! /bin/bash

Lx=1.50
nx=40
np=8
ndim=3
D_heat_g=0.1
cp=1.00
alpha=0.125
sph_c_g=10.0
prob=1.0
time_k=1.00

dname=data-nx${nx}-ndim${ndim}-Lx${Lx}-np${np}-D_heat_g${D_heat_g}-alpha${alpha}-sph_c_g${sph_c_g}-cp${cp}-prob${prob}-time_k${time_k}grad
rm -rf ${dname}
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp \
    -var alpha ${alpha} -var D_heat_g ${D_heat_g} -var ndim ${ndim} -var nx ${nx} \
    -var sph_c_g ${sph_c_g} -var cp ${cp} -var prob ${prob} -var time_k ${time_k} \
    -var Lx ${Lx} \
    -var dname ${dname}

