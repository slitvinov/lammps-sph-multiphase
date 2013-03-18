#! /bin/bash

set -e
set -u
Lx=1.0
nx=40
ndim=2
np=8
D_heat_v=0.6
D_heat_l=0.2
sph_eta_v=0.69395
cv_v=0.04
cv_l=0.04
dT=$1
Hwv=8.0
alpha=500
dprob=0.01
sph_rho_v=0.01
time_k=1.0
# parameters for kana
lmp=../../../../src/lmp_linux
mpirun=mpirun
proc="-np ${np}"

dname=data-nx${nx}-ndim${ndim}-Lx${Lx}-D_heat_v${D_heat_v}-alpha${alpha}\
-Hwv${Hwv}-dprob${dprob}-time_k${time_k}-cv_v${cv_v}-sph_rho_v${sph_rho_v}-dT${dT}\
-cv_l${cv_l}-D_heat_l${D_heat_l}-sph_eta_v${sph_eta_v}m

rm -rf ${dname}
mkdir -p ${dname}
${mpirun} ${proc} ${lmp} -in insert.lmp \
    -var alpha ${alpha} -var D_heat_v ${D_heat_v} -var ndim ${ndim} -var nx ${nx} \
    -var Hwv ${Hwv} -var dprob ${dprob} -var time_k ${time_k} \
    -var Lx ${Lx} -var cv_v ${cv_v} -var sph_rho_v ${sph_rho_v} -var dT ${dT} \
    -var cv_l ${cv_l} -var D_heat_l ${D_heat_l} \
    -var dname ${dname} -var sph_eta_v ${sph_eta_v}
