#! /bin/bash

dname=data
lmp=../../../../src/lmp_mpi
${lmp} -in poiseuille.lmp -var dname ${dname}

