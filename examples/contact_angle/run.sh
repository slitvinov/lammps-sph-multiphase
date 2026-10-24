#! /bin/bash

nx=121
dname=data-wall-nx${nx}

# number of a case
icase=2
mkdir -p ${dname}
lmp=${LMP:-lmp_mpi}
${lmp} -in droplet.lmp -var icase ${icase} -var dname ${dname} -var nx ${nx}
