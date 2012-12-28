#! /bin/bash

set -e
set -u
Lx=1.0
nx=30
ndim=2
np=1
D_heat_d=3.0
D_heat_g=1.0
sph_eta_d=0.69395
cv_d=1.0
cv_g=4.0
dT=0.1
Hwv=1.0
alpha=100
<<<<<<< HEAD
dprob=0.01
=======
dprob=0.1
>>>>>>> fae8f7622deccdd3abf64a5844ccc92d99db2ada
sph_rho_d=0.01
time_k=1.0
# parameters for kana
lmp=../../../../src/lmp_linux
mpirun=mpirun
proc="-np ${np}"

dname=data-nx${nx}-ndim${ndim}-Lx${Lx}-D_heat_d${D_heat_d}-alpha${alpha}\
-Hwv${Hwv}-dprob${dprob}-time_k${time_k}-cv_d${cv_d}-sph_rho_d${sph_rho_d}-dT${dT}\
-cv_g${cv_g}-D_heat_g${D_heat_g}-sph_eta_d${sph_eta_d}bug

rm -rf ${dname}
mkdir -p ${dname}
${lmp} -in insert.lmp \
    -var alpha ${alpha} -var D_heat_d ${D_heat_d} -var ndim ${ndim} -var nx ${nx} \
    -var Hwv ${Hwv} -var dprob ${dprob} -var time_k ${time_k} \
    -var Lx ${Lx} -var cv_d ${cv_d} -var sph_rho_d ${sph_rho_d} -var dT ${dT} \
    -var cv_g ${cv_g} -var D_heat_g ${D_heat_g} \
    -var dname ${dname} -var sph_eta_d ${sph_eta_d}
