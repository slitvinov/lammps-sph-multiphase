#! /bin/bash

set -e
set -u

nx=30
lmp=${LMP:-lmp_mpi}
${lmp} -var nx ${nx} -var dname data-nx${nx} -in bubble.lmp
