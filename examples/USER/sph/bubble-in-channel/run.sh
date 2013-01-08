#! /bin/bash

set -e
set -u
Lx=8.0
nx=20
ndim=3
np=8
D_heat_d=1.0
D_heat_g=0.1
sph_eta_d=0.69395
cv_d=0.4
cv_g=0.04
dT=0.0
Hwv=20.0
alpha=100
dprob=1.0
sph_rho_d=0.01
time_k=1.0
gy=1e2
# parameters for kana
lmp=../../../../src/lmp_linux
mpirun=mpirun
proc="-np ${np}"

dname=data-nx${nx}-ndim${ndim}-Lx${Lx}-D_heat_d${D_heat_d}-alpha${alpha}\
-Hwv${Hwv}-dprob${dprob}-time_k${time_k}-cv_d${cv_d}-sph_rho_d${sph_rho_d}-dT${dT}\
-cv_g${cv_g}-D_heat_g${D_heat_g}-sph_eta_d${sph_eta_d}-gy${gy}

rm -rf ${dname}
mkdir -p ${dname}
${mpirun} ${proc} ${lmp} -in insert.lmp \
    -var alpha ${alpha} -var D_heat_d ${D_heat_d} -var ndim ${ndim} -var nx ${nx} \
    -var Hwv ${Hwv} -var dprob ${dprob} -var time_k ${time_k} \
    -var Lx ${Lx} -var cv_d ${cv_d} -var sph_rho_d ${sph_rho_d} -var dT ${dT} \
    -var cv_g ${cv_g} -var D_heat_g ${D_heat_g} \
    -var dname ${dname} -var sph_eta_d ${sph_eta_d} -var gy ${gy}
