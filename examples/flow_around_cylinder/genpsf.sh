#! /bin/bash

set -e
set -u

datafile="$1"/flow.dat
output=${datafile/.dat/.psf}

tmpfile=$(mktemp /tmp/XXXXX)
# generate psf file with vmd
vmd -dispdev text -eofexit <<EOF
package require topotools
topo readlammpsdata  ${datafile} full
animate write psf ${tmpfile}
EOF

sed 's/ NAMD EXT//g' ${tmpfile} > ${output}
