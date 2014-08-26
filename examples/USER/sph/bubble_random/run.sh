#! /bin/bash

set -e
set -u

ndim=3
lmp=../../../../src/lmp_linux

dname=data-ndim${ndim}
/scratch/prefix-mpich3.1.2/bin/mpirun  -n 8 ${lmp} -var ndim ${ndim} -var dname ${dname} -in bubble.lmp
    
    
