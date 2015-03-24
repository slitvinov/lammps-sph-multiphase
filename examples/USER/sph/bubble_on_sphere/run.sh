#! /bin/bash

set -e
set -u

lmp=../../../../src/lmp_mpi
mpirun=mpirun
ndim=2
dname=data

${lmp} -in bubble.lmp -var ndim ${ndim} -var dname ${dname} 
