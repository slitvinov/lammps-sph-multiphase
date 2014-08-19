#! /bin/bash

rm last*
n=0
for d in $(ls -1td data*) ; do
    let "n++"
    ln -s $d last$n
done
