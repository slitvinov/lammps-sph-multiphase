#! /bin/bash

nx=121
dname=data-wall-nx${nx}

# number of a case
icase=2
mkdir -p ${dname}
lmp=../../../../src/lmp_linux
${lmp} -in droplet.lmp -var icase ${icase} -var dname ${dname} -var nx ${nx}
