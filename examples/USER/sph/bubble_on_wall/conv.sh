#! /bin/bash

for tga in *.tga; do
    printf "converting: %s\n" $tga
    convert $tga ${tga/.tga/.png}
done

convert *.png dep.gif
convert dep.gif -trim  m.gif

