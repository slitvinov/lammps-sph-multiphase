#! /bin/bash

set -e
set -u
# transform lammps dump files into punto format
awk '(NR>1)&&(FNR==1){printf "\n"; f=0} f; /ITEM: ATOMS/{f=1}' $* 
