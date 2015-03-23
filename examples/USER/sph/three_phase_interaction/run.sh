#! /bin/bash

nx=62
dname=data-wall-nx${nx}

# case number
icase=3

mkdir -p ${dname}
lmp=../../../../src/lmp_linux
${lmp} -in initial.lmp -var icase ${icase} -var dname ${dname} -var nx ${nx}

awk -f build_initial.awk ${dname}/initial.dat > ${dname}/initial.tmp
mv ${dname}/initial.tmp ${dname}/initial.dat

${lmp} -in droplet.lmp -var icase ${icase} -var dname ${dname} -var nx ${nx}
