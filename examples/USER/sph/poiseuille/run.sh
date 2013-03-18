#! /bin/bash

dname=data-wall
mpirun -np 4  ../../../../src/lmp_linux -in poiseuille.lmp -var dname ${dname}

