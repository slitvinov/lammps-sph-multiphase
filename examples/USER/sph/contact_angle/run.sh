#! /bin/bash

nx=56
dname=data-wall-nx${nx}

# case number
icase=2
mkdir -p ${dname}
lmp=../../../../src/lmp_mpi
${lmp} -in droplet.lmp -var icase ${icase} -var dname ${dname} -var nx ${nx}
