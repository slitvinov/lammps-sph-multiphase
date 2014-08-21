#! /bin/bash

nx=20
ndim=2
dname=data-nx${nx}-ndim${ndim}
lmp=../../../../src/lmp_linux
${lmp} -in insert.lmp -var ndim ${ndim} -var nx ${nx} -var dname ${dname}
