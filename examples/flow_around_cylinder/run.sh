#! /bin/bash

dname=data-wall
lmp=${LMP:-lmp_mpi}
mkdir -p ${dname}

${lmp} -in flow.lmp -var dname ${dname}

