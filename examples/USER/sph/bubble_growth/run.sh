#! /bin/bash

set -e
set -u

np=8
ndim=2
lmp=../../../../src/lmp_linux
mpirun=mpirun
proc="-np ${np}"

dname=data-ndim${ndim}

rm -rf ${dname}
mkdir -p ${dname}
${mpirun} ${proc} ${lmp} -var ndim ${ndim} -var dname ${dname} -in bubble.lmp
    
    
