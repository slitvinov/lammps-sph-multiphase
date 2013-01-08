#! /bin/bash
# filter xyz file in derectory $1

awk -f filterxyz.awk $1/data.xyz > $1/bubble.xyz

