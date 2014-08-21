#! /bin/bash

dname=data-wall
lmp=../../../../src/lmp_linux
mkdir -p ${dname}

${lmp} -in flow.lmp -var dname ${dname}

