#! /bin/bash

if [ $# -eq 2 ]; then
    n=$2
else
    n=1
fi

dn=$(ls -1td data* | awk -v n=$n 'NR==n' )
bash post.sh ${dn}

if [ $1 == "b" ]; then
    punto -B 0:0:0:1:1:1 -G 0.95:1.00 -z 1:2:3:7 ${dn}/bubble.dat
elif [ $1 == "s" ]; then
    punto -B 0:0:0:1:1:1 -G 0.95:1.00 -z 1:2:3:7 ${dn}/slice.x.dat
else
    punto -B 0:0:0:1:1:1 -G 0.95:1.00 -z 1:2:3:7 ${dn}/temp.dat
fi

echo ${dn}
