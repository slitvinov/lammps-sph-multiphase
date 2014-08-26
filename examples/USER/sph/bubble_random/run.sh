#! /bin/bash

set -e
set -u

ndim=2
lmp=../../../../src/lmp_linux

dname=data-ndim${ndim}
~/prefix-mpich/bin/mpirun -n 2 ${lmp} -var ndim ${ndim} -var dname ${dname} -in bubble.lmp
    
    
