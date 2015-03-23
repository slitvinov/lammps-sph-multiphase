#! /bin/bash

nx=40
ndim=2
dname=data-nx${nx}-ndim${ndim}
lmp=../../../../src/lmp_linux
${lmp} -in droplet.lmp -var ndim ${ndim} -var nx ${nx} -var dname ${dname}
