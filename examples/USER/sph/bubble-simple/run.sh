#! /bin/bash


Lx=1.0
nx=40
ndim=3
np=8
D_heat_g=0.05
cp=0.1
alpha=0.5
prob=0.1
sph_rho_d=1.0
time_k=1.00
cv_g=2.0

dname=data-nx${nx}-ndim${ndim}-Lx${Lx}-D_heat_g${D_heat_g}-alpha${alpha}-cp${cp}-prob${prob}-time_k${time_k}-cv_g${cv_g}-sph_rho_d${sph_rho_d}
rm -rf ${dname}
mkdir -p ${dname}
mpirun -np ${np}  ../../../../src/lmp_linux -in insert.lmp \
    -var alpha ${alpha} -var D_heat_g ${D_heat_g} -var ndim ${ndim} -var nx ${nx} \
    -var cp ${cp} -var prob ${prob} -var time_k ${time_k} \
    -var Lx ${Lx} -var cv_g ${cv_g} -var sph_rho_d ${sph_rho_d} \
    -var dname ${dname}
