#! /bin/bash

set -e
set -u

tmpfile=$(mktemp /tmp/XXXXX)
output=poly3d.psf
# generate psf file with vmd
vmd -dispdev text -eofexit <<EOF
package require topotools
topo readlammpsdata  poly3.txt full
animate write psf ${tmpfile}
EOF

sed 's/ NAMD EXT//g' ${tmpfile} > ${output}
