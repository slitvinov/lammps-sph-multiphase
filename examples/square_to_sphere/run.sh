#! /bin/bash

nx=40
ndim=3
dname=data-nx${nx}-ndim${ndim}
lmp=${LMP:-lmp_mpi}
${lmp} -in droplet.lmp -var ndim ${ndim} -var nx ${nx} -var dname ${dname}
