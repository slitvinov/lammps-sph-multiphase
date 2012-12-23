#! /bin/bash

set -e
set -u
Lx=1.0
nx=$1
ndim=2
np=8
D_heat_g=0.04
dT=0.01
Tc=0.85
Hwv=4.0
alpha=0.125
dprob=0.1
sph_rho_g=10.0
time_k=1.00
cv_g=4.0
gy=$2
# parameters for kana
lmp=../../../../src/lmp_linux
mpirun=mpirun
proc="-np ${np}"

dname=data-nx${nx}-ndim${ndim}-Lx${Lx}-D_heat_g${D_heat_g}-alpha${alpha}\
-Hwv${Hwv}-dprob${dprob}-time_k${time_k}-cv_g${cv_g}-sph_rho_g${sph_rho_g}-dT${dT}\
-gy${gy}-Tc${Tc}dr

rm -rf ${dname}
mkdir -p ${dname}
${mpirun} ${proc} ${lmp} -in insert.lmp \
    -var alpha ${alpha} -var D_heat_g ${D_heat_g} -var ndim ${ndim} -var nx ${nx} \
    -var Hwv ${Hwv} -var dprob ${dprob} -var time_k ${time_k} \
    -var Lx ${Lx} -var cv_g ${cv_g} -var sph_rho_g ${sph_rho_g} -var dT ${dT} \
    -var gy ${gy} -var Tc ${Tc} \
    -var dname ${dname}
