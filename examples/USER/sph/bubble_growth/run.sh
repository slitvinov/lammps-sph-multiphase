#! /bin/bash

set -e
set -u

np=8
ndim=3
lmp=../../../../src/lmp_linux
mpirun=mpirun
proc="-np ${np}"

dname=data-ndim${ndim}
${mpirun} ${proc} ${lmp} -var ndim ${ndim} -var dname ${dname} -in bubble.lmp
    
    
