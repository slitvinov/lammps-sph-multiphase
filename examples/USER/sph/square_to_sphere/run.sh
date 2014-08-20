#! /bin/bash

nx=20
ndim=3
dname=data-nx${nx}-ndim${ndim}
../../../../src/lmp_linux -in insert.lmp -var ndim ${ndim} -var nx ${nx} -var dname ${dname}
