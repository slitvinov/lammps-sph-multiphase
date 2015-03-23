#! /bin/bash

dname=data
lmp=../../../../src/lmp_linux
${lmp} -in poiseuille.lmp -var dname ${dname}

