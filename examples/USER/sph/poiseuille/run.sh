#! /bin/bash

dname=data
mpirun -np 4  ../../../../src/lmp_linux -in poiseuille.lmp -var dname ${dname}

