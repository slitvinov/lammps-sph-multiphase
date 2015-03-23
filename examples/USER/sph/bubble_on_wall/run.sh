#! /bin/bash

set -e
set -u

lmp=../../../../src/lmp_mpi
mpirun=mpirun
dname=data
ndim=2

${lmp} -in bubble.lmp -var ndim ${ndim} -var dname ${dname} 
