#! /bin/bash

dname=data
lmp=${LMP:-lmp_mpi}
${lmp} -in poiseuille.lmp -var dname ${dname}

