#! /bin/bash

Tc=0.0
#B="0:0:0:.75:3:2"
B="0:0:0:2:2:2"
if [ $# -eq 2 ]; then
    n=$2
else
    n=1
fi

dn=$(ls -1td data* | awk -v n=$n 'NR==n' )
bash post.sh ${dn}

if [ $1 == "a" ]; then
    cd ${dn}
    punto -B ${B} -G ${Tc}:1.00 -s 2 -z 1:2:3:7 temp.dat
elif [ $1 == "b" ]; then
    cd ${dn}
    punto -B ${B} -G ${Tc}:1.00 -z 1:2:3:7 bubble.dat
elif [ $1 == "s" ]; then
    cd ${dn}
    punto -B ${B} -G ${Tc}:1.00 -z 1:2:3:7 slice.x.dat
elif [ $1 == "c" ]; then
    cd ${dn}
    punto -B ${B} -c -z 1:2:3:8 -s 4 temp.dat
fi

echo ${dn}
