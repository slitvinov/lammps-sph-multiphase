#! /bin/bash

Tc=0.85
B="0:0:0:1.5:1.5:1.5"
if [ $# -eq 2 ]; then
    n=$2
else
    n=1
fi

dn=$(ls -1td data* | awk -v n=$n 'NR==n' )
bash post.sh ${dn}

if [ $1 == "a" ]; then
    punto -B ${B} -G ${Tc}:1.00 -s 1 -z 1:2:3:7 ${dn}/temp.dat
elif [ $1 == "b" ]; then
    punto -B ${B} -G ${Tc}:1.00 -z 1:2:3:7 ${dn}/bubble.dat
elif [ $1 == "s" ]; then
    punto -B ${B} -G ${Tc}:1.00 -z 1:2:3:7 ${dn}/slice.x.dat
elif [ $1 == "c" ]; then
    punto -B ${B} -c -z 1:2:3:8 -s 4 ${dn}/temp.dat
fi

echo ${dn}
