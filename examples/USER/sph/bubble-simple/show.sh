#! /bin/bash

Tc=0.85
if [ $# -eq 2 ]; then
    n=$2
else
    n=1
fi

dn=$(ls -1td data* | awk -v n=$n 'NR==n' )
bash post.sh ${dn}

if [ $# -eq 0 ]; then
    punto -B 0:0:0:1:1:1 -G ${Tc}:1.00 -z 1:2:3:7 ${dn}/temp.dat
elif [ $1 == "b" ]; then
    punto -B 0:0:0:1:1:1 -G ${Tc}:1.00 -z 1:2:3:7 ${dn}/bubble.dat
elif [ $1 == "s" ]; then
    punto -B 0:0:0:1:1:1 -G ${Tc}:1.00 -z 1:2:3:7 ${dn}/slice.x.dat
fi

echo ${dn}
