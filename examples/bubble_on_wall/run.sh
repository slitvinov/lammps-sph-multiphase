#! /bin/bash

set -e
set -u

lmp=${LMP:-lmp_mpi}
mpirun=mpirun
dname=data

${lmp} -in bubble.lmp -var dname ${dname} 
