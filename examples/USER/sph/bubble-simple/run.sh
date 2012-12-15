#! /bin/bash


Lx=1.5
nx=40
ndim=3
np=8
D_heat_g=0.05
cp=0.05
alpha=0.5
prob=1.0
time_k=1.00
cv_g=1.0

dname=data-nx${nx}-ndim${ndim}-Lx${Lx}-D_heat_g${D_heat_g}-alpha${alpha}-cp${cp}-prob${prob}-time_k${time_k}-cv_g${cv_g}
rm -rf ${dname}
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp \
    -var alpha ${alpha} -var D_heat_g ${D_heat_g} -var ndim ${ndim} -var nx ${nx} \
    -var cp ${cp} -var prob ${prob} -var time_k ${time_k} \
    -var Lx ${Lx} -var cv_g ${cv_g} \
    -var dname ${dname}
