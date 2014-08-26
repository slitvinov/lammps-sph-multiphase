#! /bin/bash

set -e
set -u

ndim=3
nx=20
lmp=../../../../src/lmp_linux

dname=data-ndim${ndim}-nx${nx}
${lmp} -var nx ${nx} -var ndim ${ndim} -var dname ${dname} -in bubble.lmp
    
    
