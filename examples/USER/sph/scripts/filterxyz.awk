#!/usr/bin/awk -f
# filter atoms and xyz file
BEGIN {
    # keep only B type
    pat="/[BC]/"
}

function output(    i) {
    print n
    print atm
    print atm > "/dev/stderr"
    for (i=0; i<n; i++) {
	print aux[i]
    }
    n=0
}

NF==1 && NR>1 {
    output()
    next
}

/^Atoms/{
    atm=$0
}

($1~"B")||($1~"C"){
    aux[n++]=$0
}

END {
    output()
}
