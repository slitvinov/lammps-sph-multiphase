#! /bin/bash

set -e
set -u
dTarray[160]=0.0141
dTarray[80]=0.0200
dTarray[40]=0.0283
dTarray[20]=0.0400
Lx=1.0
nx=40
ndim=3
np=8
D_heat_g=0.04
dT=0.01
cp=1.0
alpha=0.25
dprob=0.1
sph_rho_d=1.0
time_k=1.00
cv_g=4.0
# parameters for kana
lmp=../../../../src/lmp_linux
mpirun=mpirun
proc="-np ${np}"

dname=data-nx${nx}-ndim${ndim}-Lx${Lx}-D_heat_g${D_heat_g}-alpha${alpha}\
-cp${cp}-dprob${dprob}-time_k${time_k}-cv_g${cv_g}-sph_rho_d${sph_rho_d}-dT${dT}heat

rm -rf ${dname}
mkdir -p ${dname}
${mpirun} ${proc} ${lmp} -in insert.lmp \
    -var alpha ${alpha} -var D_heat_g ${D_heat_g} -var ndim ${ndim} -var nx ${nx} \
    -var cp ${cp} -var dprob ${dprob} -var time_k ${time_k} \
    -var Lx ${Lx} -var cv_g ${cv_g} -var sph_rho_d ${sph_rho_d} -var dT ${dT} \
    -var dname ${dname}
