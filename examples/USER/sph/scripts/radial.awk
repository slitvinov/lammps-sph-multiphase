#!/usr/bin/awk -f

# Make a histogram in radial direction relative to the center of the
# simulation domain

BEGIN {
    # a step of the histogram
    Rstep = 0.01
    # min radius
    Rmin = 0.0
}

inatom && $(tidx)==1 {
    r=sqrt(($(xidx)-xc)^2+($(yidx)-yc)^2+($(zidx)-zc)^2)
    histidx = int((r - Rmin)/Rstep)
    
    Tsum[histidx]+= $(Tidx)
    np[histidx]++
}

!inatom && /^ITEM: BOX / {
    getline
    xlo=$1; xhi=$2
    getline
    ylo=$1; yhi=$2
    getline
    zlo=$1; zhi=$2
    
    xc=0.5*(xlo+xhi)
    yc=0.5*(ylo+yhi)
    zc=0.5*(zlo+zhi)
}

!inatom && /^ITEM: ATOMS/ {
    # fill array with indexes
    for (i=3; i<=NF; i++) {
	idx[$(i)]=i-2
    }
    xidx=idx["x"]; 
    yidx=idx["y"]
    zidx=idx["z"]
    Tidx=idx["c_it_atom"]
    tidx=idx["type"]
    
    inatom=1
}

END {
    for (histidx in np) {
	print (histidx+1/2)*Rstep + Rmin, Tsum[histidx]/np[histidx] | "/usr/bin/sort -g"
    }
}
