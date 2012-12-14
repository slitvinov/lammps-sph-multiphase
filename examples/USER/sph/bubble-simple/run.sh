#! /bin/bash

Lx=1.0
nx=10
np=1
ndim=2
D_heat_g=0.005
cp=0.2
alpha=0.0001
sph_c_g=10.0
prob=0.0
time_k=1.8

dname=data-nx${nx}-ndim${ndim}-Lx${Lx}-np${np}-D_heat_g${D_heat_g}-alpha${alpha}-sph_c_g${sph_c_g}-cp${cp}-prob${prob}-time_k${time_k}dist
rm -rf ${dname}
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp \
    -var alpha ${alpha} -var D_heat_g ${D_heat_g} -var ndim ${ndim} -var nx ${nx} \
    -var sph_c_g ${sph_c_g} -var cp ${cp} -var prob ${prob} -var time_k ${time_k} \
    -var Lx ${Lx} \
    -var dname ${dname}

