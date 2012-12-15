#! /bin/bash

dn=$(ls -1trd data* | tail -1)
bash post.sh ${dn}

if [ $1 == "b" ]; then
    punto -B 0:0:0:1:1:1 -G 0.95:1.00 -z 1:2:3:7 ${dn}/bubble.dat
elif [ $1 == "s" ]; then
    punto -B 0:0:0:1:1:1 -G 0.95:1.00 -z 1:2:3:7 ${dn}/slice.x.dat
elif [ $1 == "c" ]; then
    punto -B 0:0:0:1:1:1 -c -z 1:2:3:8 ${dn}/temp.dat
else
    punto -B 0:0:0:1:1:1 -G 0.95:1.00 -z 1:2:3:7 ${dn}/temp.dat
fi

echo ${dn}
