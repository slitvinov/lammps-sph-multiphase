#! /bin/bash

Lx=1.50
nx=40
np=8
ndim=2
D_heat_g=0.01
cp=10.0
alpha=0.125
prob=0.05
time_k=1.00
cv_g=10.0

dname=data-nx${nx}-ndim${ndim}-Lx${Lx}-D_heat_g${D_heat_g}-alpha${alpha}-cp${cp}-prob${prob}-time_k${time_k}-cv_g${cv_g}
rm -rf ${dname}
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp \
    -var alpha ${alpha} -var D_heat_g ${D_heat_g} -var ndim ${ndim} -var nx ${nx} \
    -var cp ${cp} -var prob ${prob} -var time_k ${time_k} \
    -var Lx ${Lx} -var cv_g ${cv_g} \
    -var dname ${dname}

