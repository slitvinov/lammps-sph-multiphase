#! /bin/bash

set -e
set -u

lmp=../../../../src/lmp_linux
mpirun=mpirun
dname=data

${lmp} -in bubble.lmp -var dname ${dname} 
