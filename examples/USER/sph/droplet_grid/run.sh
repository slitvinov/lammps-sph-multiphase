#! /bin/bash

# 40, 68, 126
nx=68
dname=data-wall-nx${nx}

# case number (1 or 2)
icase=1
mkdir -p ${dname}
../../../../src/lmp_linux -in droplet.lmp -var icase ${icase} -var dname ${dname} -var nx ${nx}
