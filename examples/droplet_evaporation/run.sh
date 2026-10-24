#! /bin/bash

set -e
set -u

ndim=3
nx=30
lmp=${LMP:-lmp_mpi}
${lmp} -var ndim ${ndim} -var nx ${nx} -var dname data-ndim${ndim}-nx${nx} -in droplet.lmp
