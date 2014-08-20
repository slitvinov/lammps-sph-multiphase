#! /bin/bash

set -e
set -u

input=$1/droplet.restart
datafile=${input/.restart/.dat}
output=${input/.restart/.psf}

../../../../tools/restart2data ${input} ${datafile}

tmpfile=$(mktemp /tmp/XXXXX)
# generate psf file with vmd
vmd -dispdev text -eofexit <<EOF
package require topotools
topo readlammpsdata  ${datafile} full
animate write psf ${tmpfile}
EOF

sed 's/ NAMD EXT//g' ${tmpfile} > ${output}
