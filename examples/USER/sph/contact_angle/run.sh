#! /bin/bash

nx=248
dname=data-wall-nx${nx}

# number of a case
icase=1
mkdir -p ${dname}
lmp=../../../../src/lmp_linux
${lmp} -in droplet.lmp -var icase ${icase} -var dname ${dname} -var nx ${nx}
